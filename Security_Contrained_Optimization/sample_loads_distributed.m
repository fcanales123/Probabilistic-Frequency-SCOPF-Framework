function [mpc_s, PD_s_MW, info] = sample_loads_distributed(mpc0, ldc, sampler)
%SAMPLE_LOADS_DISTRIBUTED Sample distributed bus loads using:
%   PD_i = m * PD0_i * (1 + eps_i), eps_i ~ N(0, sigma_i)
%
% Inputs
%   mpc0.bus(:,PD) in MW (MATPOWER convention)
%   ldc.P_MW : vector of system total load levels (MW), e.g. from LDC
%   ldc.u    : cumulative probability grid in [0,1] same length as ldc.P_MW
%   sampler fields (recommended):
%     .method        'inv_cdf' or 'hour_pick'      (for global factor m)
%     .sigma_rel     scalar or nbx1 vector of std-dev on per-bus multiplier
%                    (e.g. 0.05 means 5% std dev)
%     .clip_mult     [] or [lo hi] for per-bus multiplier (after noise), e.g. [0.5 1.5]
%     .clip_m        [] or [lo hi] for global multiplier m (optional)
%     .keep_sum      true/false renormalize PD to match m*sum(PD0)
%     .scale_QD      true/false scale QD proportionally to PD
%     .seed          [] or integer seed for reproducibility
%
% Outputs
%   mpc_s       : case with sampled bus PD (and optionally QD) in MW
%   PD_s_MW     : nbx1 sampled real demands in MW
%   info        : struct with diagnostics

    define_constants;

    if ~isfield(sampler, 'method'),    sampler.method = 'inv_cdf'; end
    if ~isfield(sampler, 'sigma_rel'), sampler.sigma_rel = 0.05; end
    if ~isfield(sampler, 'clip_mult'), sampler.clip_mult = []; end
    if ~isfield(sampler, 'clip_m'),    sampler.clip_m = []; end
    if ~isfield(sampler, 'keep_sum'),  sampler.keep_sum = true; end
    if ~isfield(sampler, 'scale_QD'),  sampler.scale_QD = true; end
    if ~isfield(sampler, 'seed'),      sampler.seed = []; end

    if ~isempty(sampler.seed)
        rng(sampler.seed);
    end

    % --- base loads (MW) ---
    PD0_MW = mpc0.bus(:, PD);
    QD0_MW = mpc0.bus(:, QD);
    nb = size(mpc0.bus, 1);

    Psys0_MW = sum(PD0_MW);

    % --- LDC support vectors ---
    P = ldc.P_MW(:);
    ugrid = ldc.u(:);

    if numel(P) ~= numel(ugrid)
        error('sample_loads_distributed: ldc.P_MW and ldc.u must have same length.');
    end

    % --- sample system total load from LDC ---
    switch lower(sampler.method)
        case 'inv_cdf'
            u = rand();
            Psys_MW = interp1(ugrid, P, u, 'linear', 'extrap');
        case 'hour_pick'
            idx = randi(numel(P));
            Psys_MW = P(idx);
        otherwise
            error('sample_loads_distributed: unknown sampler.method "%s".', sampler.method);
    end

    if ~isfinite(Psys_MW) || Psys_MW < 0
        error('sample_loads_distributed: sampled Psys_MW invalid (%.4g).', Psys_MW);
    end

    % --- global multiplier m ---
    m = Psys_MW / Psys0_MW;

    % optional clipping on m
    if ~isempty(sampler.clip_m)
        m = min(max(m, sampler.clip_m(1)), sampler.clip_m(2));
        Psys_MW = m * Psys0_MW;
    end

    % --- per-bus noise std-dev(s) ---
    if isscalar(sampler.sigma_rel)
        sigma = sampler.sigma_rel * ones(nb, 1);
    else
        sigma = sampler.sigma_rel(:);
        if numel(sigma) ~= nb
            error('sample_loads_distributed: sigma_rel must be scalar or nbx1.');
        end
    end

    % Only apply noise to buses with positive base load
    active = PD0_MW > 0;

    % Gaussian noise on multipliers: mult_i = 1 + N(0, sigma_i)
    eps_i = zeros(nb, 1);
    eps_i(active) = sigma(active) .* randn(nnz(active), 1);

    mult_i = ones(nb, 1);
    mult_i(active) = 1 + eps_i(active);

    % optional clipping on per-bus multiplier
    if ~isempty(sampler.clip_mult)
        lo = sampler.clip_mult(1);
        hi = sampler.clip_mult(2);
        mult_i(active) = min(max(mult_i(active), lo), hi);
    end

    % preliminary sampled PD (MW)
    PD_s_MW = zeros(nb, 1);
    PD_s_MW(active) = m * PD0_MW(active) .* mult_i(active);

    % enforce nonnegativity
    PD_s_MW(PD_s_MW < 0) = 0;

    % --- optional renormalization to hit target system total exactly ---
    if sampler.keep_sum
        target = m * sum(PD0_MW);      % = Psys_MW (after clip_m)
        cur = sum(PD_s_MW);
        if cur > 0
            scale = target / cur;
            PD_s_MW(active) = PD_s_MW(active) * scale;
        else
            scale = NaN;
        end
    else
        scale = NaN;
    end

    % --- build output case ---
    mpc_s = mpc0;
    mpc_s.bus(:, PD) = PD_s_MW;

    if sampler.scale_QD
        % scale QD proportional to PD at each bus (keep power factor)
        ratio = ones(nb, 1);
        ratio(active) = PD_s_MW(active) ./ PD0_MW(active);
        mpc_s.bus(:, QD) = QD0_MW .* ratio;
    end

    % --- outputs / diagnostics ---
    info = struct();
    info.Psys0_MW = Psys0_MW;
    info.Psys_MW  = sum(PD_s_MW);
    info.Psys_draw_MW = Psys_MW;
    info.m_global = m;
    info.mult_i_mean = mean(mult_i(active));
    info.mult_i_std  = std(mult_i(active));
    info.renorm_scale = scale;
    info.method = sampler.method;
end
