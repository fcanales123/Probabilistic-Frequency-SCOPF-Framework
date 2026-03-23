function [mpc_s, mult, Psys_MW, info] = sample_loads(mpc0, ldc, sampler)
%SAMPLE_LOADS  Sample system load from an LDC and scale MATPOWER bus loads.
%
%   [mpc_s, mult, Psys_MW, info] = sample_loads(mpc0, ldc, sampler)
%
% Inputs
%   mpc0    : base MATPOWER case (external indexing OK)
%   ldc     : struct with fields:
%               .P_MW    (8760x1) system load values in MW (your scaled_load)
%               .u       (8760x1) CDF support in [0,1] (your total_hours)
%             OR optionally:
%               .hours   (8760x1) optional, not used
%   sampler : struct (optional) controls behavior:
%               .method  = 'inv_cdf' (default) or 'hour_pick'
%               .rng     = []  (optional) set rng seed externally if you want
%               .clip    = [0.5 1.5] (optional) clamp multiplier
%               .keepPd0 = true (default) store PD0 in info
%
% Outputs
%   mpc_s   : sampled MATPOWER case with bus(:,PD) scaled by mult
%   mult    : system-wide multiplier applied to all PD
%   Psys_MW : sampled system total load in MW
%   info    : struct with diagnostics

    define_constants;

    if nargin < 2 || isempty(ldc)
        error('sample_loads: ldc input is required with fields .P_MW and .u');
    end
    if nargin < 3, sampler = struct(); end

    if ~isfield(sampler, 'method') || isempty(sampler.method)
        sampler.method = 'inv_cdf';
    end
    if ~isfield(sampler, 'clip'), sampler.clip = []; end
    if ~isfield(sampler, 'keepPd0'), sampler.keepPd0 = true; end

    % --- extract base PD pattern (MW) ---
    PD0_MW = mpc0.bus(:, PD);
    Psys0_MW = sum(PD0_MW);

    if Psys0_MW <= 0
        error('sample_loads: base case has non-positive total load (sum PD <= 0).');
    end

    % --- validate LDC vectors ---
    P = ldc.P_MW(:);
    if isfield(ldc, 'u')
        ugrid = ldc.u(:);
    else
        % fallback: assume uniform CDF grid
        ugrid = linspace(0, 1, numel(P)).';
    end

    if numel(P) ~= numel(ugrid)
        error('sample_loads: ldc.P_MW and ldc.u must have same length.');
    end
    if any(~isfinite(P)) || any(~isfinite(ugrid))
        error('sample_loads: ldc vectors contain NaN/Inf.');
    end

    % Ensure ugrid is strictly increasing for interp1
    % (Your total_hours is monotone increasing, so this should be fine.)
    [ugrid, iu] = unique(ugrid, 'stable');
    P = P(iu);

    % --- sample Psys_MW ---
    switch lower(sampler.method)
        case 'inv_cdf'
            % Inverse transform sampling from empirical CDF
            u = rand();
            % Interpolate on u in [0,1]; extrapolate just in case of numerical edges
            Psys_MW = interp1(ugrid, P, u, 'linear', 'extrap');

        case 'hour_pick'
            % Equivalent to sampling a random hour from the year (empirical distribution)
            idx = randi(numel(P));
            Psys_MW = P(idx);

        otherwise
            error('sample_loads: unknown sampler.method "%s".', sampler.method);
    end

    % guard
    if ~isfinite(Psys_MW) || Psys_MW < 0
        error('sample_loads: sampled Psys_MW is invalid (%.4g).', Psys_MW);
    end

    % --- compute multiplier ---
    mult = Psys_MW / Psys0_MW;

    % Optional clipping (useful to avoid insane draws while debugging)
    if ~isempty(sampler.clip)
        lo = sampler.clip(1);
        hi = sampler.clip(2);
        mult = min(max(mult, lo), hi);
        Psys_MW = mult * Psys0_MW;
    end

    % --- apply scaling ---
    mpc_s = mpc0;
    mpc_s.bus(:, PD) = mult * PD0_MW;

    % --- outputs / diagnostics ---
    info = struct();
    info.Psys0_MW = Psys0_MW;
    info.Psys_MW  = Psys_MW;
    info.mult     = mult;
    info.method   = sampler.method;
    if sampler.keepPd0
        info.PD0_MW = PD0_MW;
    end
end
