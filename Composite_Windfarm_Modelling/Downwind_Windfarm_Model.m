function [farm, perBin, diag] = Downwind_Windfarm_Model(X, Y, dir_bins_deg, speed_bins, prob, turbine, wT_curve, params)
% Downwind_Windfarm_Model
% Propagate-downwind engineering model (Bastankhah Gaussian + RSS superposition)
% to compute per-turbine inflow wind speeds and wind-farm power for each
% (direction, speed) bin, and the expectation over the joint distribution.
%
% Inputs
%   X, Y               : column vectors (nT x 1) turbine coordinates in meters
%   dir_bins_deg       : vector of wind directions [deg FROM which wind blows]
%   speed_bins         : vector of reference hub-height wind speeds [m/s]
%   prob               : probability weights. Either
%                        - matrix size [numDirs x numSpeeds] summing to 1, or
%                        - scalar 1 (uniform assumed), or empty [] (uniform)
%   turbine            : 'IEA' or 'DTU' (used for rotor diameter if not in params)
%   wT_curve           : numeric array with columns [WS, Power, CT, CP]
%                        Units: WS [m/s], Power [W] (or kW/MW — see params.power_unit),
%                        CT [-], CP [-]. If Power column is unavailable or 0,
%                        the model can compute power from CP.
%   params             : struct with optional fields
%                        .rho (kg/m^3, default 1.225)
%                        .k, .ceps, .ctlim for Bastankhah (defaults align with your fn)
%                        .D (m) rotor diameter (if omitted, inferred from turbine type)
%                        .power_unit ('W'|'kW'|'MW'), default 'W'
%                        .use_power_curve (true/false), default true
%                        .cut_in, .cut_out (m/s) optional; if absent, inferred from curve
%                        .Prated (W) optional cap if computing from CP
%
% Outputs
%   farm: struct with expected values across bins
%       .P_expected_W        : expected farm power [W]
%       .P_per_dir_W         : expected power per direction (marginalized over speeds)
%       .P_per_speed_W       : expected power per speed (marginalized over directions)
%       .notes               : modeling assumptions
%   perBin: struct with detailed arrays for each bin
%       .Pbin_W (nDir x nSpd)      : farm power per bin
%       .Uhub{d,s} (nT x 1)        : per-turbine inflow speeds per bin
%       .Pti_W{d,s} (nT x 1)       : per-turbine power per bin
%       .CTi{d,s} (nT x 1)         : per-turbine CT used per bin
%
% Dependencies: bastankhahGaussianDeficit.m available on path.
%
% Author: Fernando/ChatGPT

nT = numel(X);
X = X(:); Y = Y(:);
ndir = numel(dir_bins_deg);
ns = numel(speed_bins);

% --- Params & curves ---
if ~isfield(params, 'rho') || isempty(params.rho), params.rho = 1.225; end
if ~isfield(params, 'k'),     params.k = 0.0324555; end
if ~isfield(params, 'ceps'),  params.ceps = 0.2; end
if ~isfield(params, 'ctlim'), params.ctlim = 0.899; end
if ~isfield(params, 'power_unit') || isempty(params.power_unit), params.power_unit = 'W'; end
if ~isfield(params, 'use_power_curve') || isempty(params.use_power_curve), params.use_power_curve = true; end

% rotor diameter
if isfield(params,'D') && ~isempty(params.D)
    D = params.D;
else
    switch upper(turbine)
        case 'IEA', D = 240;
        case 'DTU', D = 178.3;
        otherwise, error('Unknown turbine type and D not provided.');
    end
end
A = pi*(D/2)^2; rho = params.rho;

% Parse wT_curve columns: [WS, Power, CT, CP]
ws_tab = wT_curve(:,1);
P_tab  = wT_curve(:,2);
CT_tab = wT_curve(:,3);
CP_tab = wT_curve(:,4);

% Unit handling for power column if used
switch lower(params.power_unit)
    case 'w'
        P_tab_W = P_tab;
    case 'kw'
        P_tab_W = P_tab * 1e3;
    case 'mw'
        P_tab_W = P_tab * 1e6;
    otherwise
        error('Unknown params.power_unit. Use ''W'',''kW'',''MW''.');
end

% Cut-in / cut-out / rated inference (optional override via params)
if isfield(params,'cut_in') && ~isempty(params.cut_in)
    cut_in = params.cut_in; else
    pos = find((CP_tab>0) | (P_tab>0), 1, 'first'); if isempty(pos), pos = 1; end
    cut_in = ws_tab(pos);
end
if isfield(params,'cut_out') && ~isempty(params.cut_out)
    cut_out = params.cut_out; else
    pos = find((CP_tab>0) | (P_tab>0), 1, 'last'); if isempty(pos), pos = numel(ws_tab); end
    cut_out = ws_tab(pos);
end
if isfield(params,'rated') && ~isempty(params.rated)
    rated = params.rated; else
    % infer rated as first occurrence of max power (robust to plateaus)
    [~, ridx] = max(P_tab);
    rated = ws_tab(ridx);
end

% Optional rated power if computing from CP
Prated = getfield_safe(params,'Prated',inf);

% --- piecewise interpolation helpers ---
% POWER: pchip for [cut_in,rated), linear for [rated,cut_out), nearest at cut_out
interp_power_piecewise = @(ws_grid, val_grid, q) ...
    local_piecewise_interp_power(ws_grid, val_grid, q, cut_in, rated, cut_out);

% CP/CT: linear for [cut_in,rated), pchip for [rated,cut_out), nearest at cut_out
interp_cpct_piecewise = @(ws_grid, val_grid, q) ...
    local_piecewise_interp_cpct(ws_grid, val_grid, q, cut_in, rated, cut_out);


% probability handling
if isempty(prob)
    prob = ones(ndir,ns)/(ndir*ns);
elseif isscalar(prob) && prob==1
    prob = ones(ndir,ns)/(ndir*ns);
else
    if ~isequal(size(prob), [ndir, ns])
        error('prob must be size [numel(dir_bins_deg) x numel(speed_bins)]');
    end
    s = sum(prob(:));
    if abs(s-1)>1e-6
        prob = prob / s; % normalize defensively
    end
end

% outputs
perBin.Pbin_W = zeros(ndir, ns);
perBin.Uhub = cell(ndir, ns);
perBin.Pti_W = cell(ndir, ns);
perBin.CTi = cell(ndir, ns);

for id = 1:ndir
    theta = dir_bins_deg(id);

    % Rotate into wind-aligned frame: x = downstream, y = crosswind
    % Wind blowing FROM theta degrees (0=x+, 90=y+). Downwind unit vector e = [cosd theta, sind theta].
    c = cosd(theta); s = sind(theta);
    xdw =  c*X + s*Y;   % projection along downwind axis
    ycw = -s*X + c*Y;   % crosswind (left-positive)

    % Order turbines from most upwind (smallest xdw) to downwind (largest xdw)
    [xdw_sorted, order] = sort(xdw, 'ascend'); %#ok<ASGLU>
    ycw_sorted = ycw(order);
    idx_inv(order) = 1:nT; idx_inv = idx_inv(:); %#ok<AGROW>

    for ispd = 1:ns
        U0 = speed_bins(ispd);

        U_i   = zeros(nT,1);    % inflow at each turbine
        P_i   = zeros(nT,1);    % power at each turbine
        CT_i  = zeros(nT,1);

        % propagate downwind
        for ii = 1:nT
            iT = order(ii);
            xi = xdw_sorted(ii);
            yi = ycw_sorted(ii);

            deficits = [];
            % Accumulate wakes from upstream machines k
            deltas = [];   % normalized deficits (1 - Uw_i/U0)
            fracs  = [];   % area overlap fractions (Aw/Ar)
            for kk = 1:ii-1
                kT = order(kk);
                xk = xdw_sorted(kk);
                yk = ycw_sorted(kk);
                dw = xi - xk;        % downstream distance (>0)
                cw = yi - yk;        % lateral offset

                Uk = U_i(kT);        % upstream inflow already computed
                if Uk<=0, Uk = 0; end

                % CT of upstream turbine at its inflow speed (from curve)
                CTk = interp_cpct_piecewise(ws_tab, CT_tab, Uk);
                CTk = max(0, min(CTk, params.ctlim));

                % Single-wake from k in undisturbed wind U0 at turbine i location
                [d_ik, wakeDiam_k, ~] = bastankhahGaussianDeficit(D, dw, cw, CTk, U0, ...
                    'k', params.k, 'ceps', params.ceps, 'ctlim', params.ctlim);
                Uwik = max(0, U0 - d_ik);
                delta_k = max(0, 1 - Uwik / U0);   % normalized deficit for mixed-wake RSS

                % Partial wake incidence: overlap of wake disk and rotor disk (uniform deficit assumption)
                Rr = D/2;                  % rotor radius
                Rw = (wakeDiam_k)/2;       % wake radius at this x-section
                Aw = circle_overlap_area(Rr, Rw, abs(cw));
                Ar = pi*Rr^2;
                frac = max(0, min(1, Aw/Ar));

                if frac>0
                    deltas(end+1,1) = delta_k; %#ok<AGROW>
                    fracs(end+1,1)  = frac;    %#ok<AGROW>
                end
            end

            % Mixed-wake combination
            if isempty(deltas)
                Ui = U0;   % no wake coverage
            else
                % Mixed-wake (within covered region): RSS of normalized deficits
                rss = sqrt(sum(deltas.^2));
                Uw_mixed = U0 * (1 - rss);

                % Area weighting for partial coverage (uniform-deficit assumption)
                Fcover = max(0, min(1, sum(fracs))); % approx. union fraction
                Ui = Fcover * Uw_mixed + (1 - Fcover) * U0;
            end
            Ui = max(0, Ui);

            % Turbine-level power: prefer provided power curve when requested
            if params.use_power_curve && any(P_tab_W>0)
                Pi = interp_power_piecewise(ws_tab, P_tab_W, Ui);
                if isnan(Pi) || Ui<cut_in || Ui>cut_out, Pi = 0; end
            else
                % compute from CP curve
                CPi = max(0, interp_cpct_piecewise(ws_tab, CP_tab, Ui));
                Pi  = 0.5 * rho * A * CPi * Ui^3;
                if ~isinf(Prated), Pi = min(Pi, Prated); end
            end

            % Apply cut-in / cut-out (zero power outside)
            if Ui<cut_in || Ui>cut_out
                Pi = 0; Ui = 0; % assume parked
            end

            % Store

            % Store
            U_i(iT)  = Ui;
            P_i(iT)  = Pi;
            CT_i(iT) = max(0, interp_cpct_piecewise(ws_tab, CT_tab, Ui));
        end

        perBin.Uhub{id,ispd} = U_i;
        perBin.Pti_W{id,ispd} = P_i;
        perBin.CTi{id,ispd} = CT_i;
        perBin.Pbin_W(id,ispd) = sum(P_i);
    end
end

% Expected values
farm.P_expected_W  = sum(perBin.Pbin_W(:) .* prob(:));
farm.P_per_dir_W   = sum(perBin.Pbin_W .* prob, 2);                  % sum over speeds
farm.P_per_speed_W = sum(perBin.Pbin_W .* prob, 1);                  % sum over dirs
farm.notes = sprintf(['Downwind model with Bastankhah Gaussian (k=%.5f, ceps=%.3f), RSS superposition. ', ...
                      'CT/P/CP from wT_curve [WS,Power,CT,CP]. Prefer power curve=%d. ', ...
                      'Cut-in=%.2f m/s, Cut-out=%.2f m/s. Directions meteorological (FROM).'], ...
                      params.k, params.ceps, params.use_power_curve, cut_in, cut_out);

% ---- Diagnostics ----
U0_grid = ones(ndir,1) * speed_bins(:)';
diag.mean_U0 = sum(U0_grid(:) .* prob(:));
% mean inflow at rotor (averaged over turbines and bins)
meanUi = 0;
for id = 1:ndir
for ispd = 1:ns
Ui_vec = perBin.Uhub{id,ispd};
if isempty(Ui_vec), continue; end
meanUi = meanUi + mean(Ui_vec) * prob(id,ispd);
end
end
diag.mean_U_inflow = meanUi;
diag.wake_speed_ratio = safe_div(meanUi, diag.mean_U0);


% Rated & specific power
if params.use_power_curve && any(P_tab_W>0)
Prated_guess = max(P_tab_W);
else
% compute from CP at rated speed if possible
CP_r = max(0, interp_piecewise(ws_tab, CP_tab, rated));
Prated_guess = 0.5 * rho * A * CP_r * rated^3;
end


diag.turbine = turbine;
diag.D = D; diag.Area_m2 = A; diag.Prated_W = Prated_guess;
diag.specific_power_Wm2 = safe_div(Prated_guess, A);


% Expected farm power with NO wakes (all turbines at U0)
P_nowake = 0;
for id = 1:ndir
for ispd = 1:ns
U0 = speed_bins(ispd);
if params.use_power_curve && any(P_tab_W>0)
Pi0 = interp_power_piecewise(ws_tab, P_tab_W, U0);
if isnan(Pi0) || U0<cut_in || U0>cut_out, Pi0 = 0; end
else
CP0 = max(0, interp_cpct_piecewise(ws_tab, CP_tab, U0));
Pi0 = 0.5 * rho * A * CP0 * U0^3;
if ~isinf(Prated), Pi0 = min(Pi0, Prated); end
if U0<cut_in || U0>cut_out, Pi0 = 0; end
end
P_nowake = P_nowake + nT * Pi0 * prob(id,ispd);
end
end


diag.ExpectedPowerNoWake_W = P_nowake;
diag.WakeLossFraction = 1 - safe_div(farm.P_expected_W, max(P_nowake, eps));
end


function r = safe_div(a,b)
if b<=0, r = NaN; else, r = a/b; end
end

function A = circle_overlap_area(R1, R2, d)
% Area of overlap of two circles of radii R1, R2 with center distance d
% Handles disjoint and full-containment cases.
if d >= R1 + R2
    A = 0; return; % no overlap
elseif d <= abs(R1 - R2)
    % one circle fully inside the other
    A = pi * min(R1,R2)^2; return;
else
    % partial overlap
    alpha = 2*acos( (d^2 + R1^2 - R2^2) / (2*d*R1) );
    beta  = 2*acos( (d^2 + R2^2 - R1^2) / (2*d*R2) );
    A = 0.5*R1^2*(alpha - sin(alpha)) + 0.5*R2^2*(beta - sin(beta));
end
end

function y = local_piecewise_interp_power(ws_grid, val_grid, q, cut_in, rated, cut_out)
% POWER: pchip [cut_in,rated), linear [rated,cut_out), nearest at cut_out
if numel(q)~=1
    y = arrayfun(@(qq) local_piecewise_interp_power(ws_grid, val_grid, qq, cut_in, rated, cut_out), q); return;
end
if q < cut_in
    y = interp1(ws_grid, val_grid, cut_in, 'nearest', 'extrap'); if y<eps, y=0; end; return;
elseif q > cut_out
    y = interp1(ws_grid, val_grid, cut_out, 'nearest', 'extrap'); if y<eps, y=0; end; return;
elseif q == cut_out
    y = interp1(ws_grid, val_grid, q, 'nearest', 'extrap'); return;
elseif q < rated
    y = interp1(ws_grid, val_grid, q, 'pchip'); return;
else
    y = interp1(ws_grid, val_grid, q, 'linear'); return;
end
end

function y = local_piecewise_interp_cpct(ws_grid, val_grid, q, cut_in, rated, cut_out)
% CP/CT: linear [cut_in,rated), pchip [rated,cut_out), nearest at cut_out
if numel(q)~=1
    y = arrayfun(@(qq) local_piecewise_interp_cpct(ws_grid, val_grid, qq, cut_in, rated, cut_out), q); return;
end
if q < cut_in
    y = interp1(ws_grid, val_grid, cut_in, 'nearest', 'extrap'); return;
elseif q > cut_out
    y = interp1(ws_grid, val_grid, cut_out, 'nearest', 'extrap'); return;
elseif q == cut_out
    y = interp1(ws_grid, val_grid, q, 'nearest', 'extrap'); return;
elseif q < rated
    y = interp1(ws_grid, val_grid, q, 'linear'); return;
else
    y = interp1(ws_grid, val_grid, q, 'pchip'); return;
end
end

function v = getfield_safe(s, f, default)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = default; end
end

%% --- Convenience runner (example) ---
%{
% Example usage:
%
% % 1) Layout
% nT = 25; turbine = 'IEA'; spacing = 7;  % 7D spacing
% [X,Y] = Generate_Honeycomb_Layout(nT, turbine, spacing);
%
% % 2) Bins & probabilities
% dir_bins = 0:30:330;              % 12 directions
% speed_bins = 4:2:24;              % 4,6,...,24 m/s
% prob = ones(numel(dir_bins), numel(speed_bins)); prob = prob/sum(prob(:));
%
% % 3) Curves loaded earlier into wT_curve (N x 4): [WS, Power, CT, CP]
% % params.power_unit = 'W'; % or 'kW' or 'MW' depending on your file
% params.D = [];        % use default from turbine type if empty
% params.rho = 1.225; params.k = 0.0324555; params.ceps = 0.2; params.ctlim = 0.899;
% params.use_power_curve = true; % set false to compute from CP
% % params.cut_in = 3; params.cut_out = 25; % optional explicit
%
% % 4) Run
% [farm, perBin] = Downwind_Windfarm_Model(X, Y, dir_bins, speed_bins, prob, turbine, wT_curve, params);
% fprintf('Expected farm power = %.2f MW
', farm.P_expected_W/1e6);
%}
