function [MC, results, CI] = mc_expected_droop_CI(mpc, ldc, sampler, mpopt, opts)
%MC_EXPECTED_DROOP_CI Monte Carlo expected droop-stage congestion index.
%
% Returns running statistics for E[ CI_droop ] where CI_droop is already
% contingency-weighted inside csopf_droop_CI via CI.E_CI = sum(p_k*CI_system(k)).

    if nargin < 4 || isempty(mpopt)
        mpopt = mpoption('model','DC','opf.dc.solver','MIPS','verbose',0);
    end
    if nargin < 5 || isempty(opts)
        opts = struct();
    end

    % --- stopping / safety defaults ---
    if ~isfield(opts,'Nmin');        opts.Nmin = 50; end
    if ~isfield(opts,'Nmax');        opts.Nmax = 5000; end
    if ~isfield(opts,'rel_tol');     opts.rel_tol = 0.01; end      % 1% relative half-width
    if ~isfield(opts,'abs_tol');     opts.abs_tol = 1e-3; end      % absolute half-width
    if ~isfield(opts,'conf');        opts.conf = 0.95; end         % CI level
    if ~isfield(opts,'print_every'); opts.print_every = 25; end
    if ~isfield(opts,'mode');        opts.mode = 'max'; end        % for CI_system scalar inside csopf_droop_CI

    define_constants;

    % --- wind sampling defaults ---
    if ~isfield(opts,'wind'); opts.wind = struct(); end
    if ~isfield(opts.wind,'enable'); opts.wind.enable = isfield(mpc,'wind') && isfield(mpc.wind,'idx'); end
    if ~isfield(opts.wind,'clip_pu'); opts.wind.clip_pu = [0 1]; end  % clip sampled pu
    if opts.wind.enable
        if ~isfield(opts.wind,'models')
            error('mc_expected_droop_CI: opts.wind.models is required when wind sampling is enabled.');
        end
        if ~isfield(mpc,'wind') || ~isfield(mpc.wind,'idx')
            error('mc_expected_droop_CI: mpc.wind.idx not found. Store wind gen indices when you add wind units.');
        end
        if ~isfield(mpc.wind,'P_rated_MW')
            error('mc_expected_droop_CI: mpc.wind.P_rated_MW not found. Store rated MW per wind unit.');
        end
    end


    % --- running stats (Welford) for E_CI across samples ---
    n = 0;
    meanE = 0;
    M2 = 0;

    % bookkeeping
    E_samples = NaN(opts.Nmax,1);
    mult_samples = NaN(opts.Nmax,1);
    success = false(opts.Nmax,1);

    have_vec = false;

    mean_RDp = [];  mean_RDm = [];
    mean_LS  = [];  mean_LFp = [];  mean_LFm = [];
    
    RD_gen_idx = [];
    LF_bus_idx = [];

    % confidence z (normal approx; fine for n>~30)
    z = 1.96;
    if abs(opts.conf - 0.90) < 1e-12, z = 1.645; end
    if abs(opts.conf - 0.95) < 1e-12, z = 1.96;  end
    if abs(opts.conf - 0.99) < 1e-12, z = 2.576; end

    for it = 1:opts.Nmax

        % ---- (1) sample WPP output and modify mpc ----
        mpc_w = mpc;   % start from base each iteration

        if opts.wind.enable
            wind_idx = mpc.wind.idx(:);
            P_rated  = mpc.wind.P_rated_MW(:);
            nw = numel(wind_idx);

            % sample pu availability (nw x 1)
            pw_pu = arrayfun(@(k) opts.wind.models(k).sample(1), (1:nw)).';

            % clip pu (KDE sometimes overshoots slightly)
            lo = opts.wind.clip_pu(1); hi = opts.wind.clip_pu(2);
            pw_pu = min(max(pw_pu, lo), hi);

            % available MW per wind unit
            P_avail_MW = pw_pu .* P_rated;

            % write into gen table
            mpc_w.gen(wind_idx, PG)   = P_avail_MW;
            mpc_w.gen(wind_idx, PMAX) = P_avail_MW;
            mpc_w.gen(wind_idx, PMIN) = 0;      % allow curtailment (optional but recommended)
        end

        if isfield(opts,'wel') && opts.wel.enable
            wel_idx   = mpc.wind_elec.idx(:);         % rows in mpc.gen of the WEL net generators
            Pw_cap    = mpc.wind_elec.Pw_cap_MW(:);   % MW capacity for each WEL's windfarm
            alpha     = mpc.wind_elec.alpha(:);
            beta      = mpc.wind_elec.beta(:);
            nwel = numel(wel_idx);
        
            % clip bounds for sampled Pw_avail_pu
            lo = opts.wel.clip_pu(1);
            hi = opts.wel.clip_pu(2);
        
            Pwel_inj_MW = zeros(nwel,1);
        
            for i = 1:nwel
                % sample wind availability PU from the correct FCOPT table for this capacity
                [Pw_avail_pu, ~, ~] = sample_wel_from_FCOPT( ...
                    opts.wel.FCOPT_cell{i}, Pw_cap(i), opts.wel.turbine{i});
        
                Pw_avail_pu = min(max(Pw_avail_pu, lo), hi);
                Pw_avail_MW = Pw_avail_pu * Pw_cap(i);
        
                % parameters
                Pel_cap_MW  = beta(i)  * Pw_cap(i);
                Pgrid_con_MW = alpha(i) * Pel_cap_MW;         % = alpha*beta*Pw_cap
        
                % operational rule
                Pel_hot_MW  = min(Pel_cap_MW, Pw_avail_MW);
                Pwel_inj_MW(i) = Pw_avail_MW - Pel_hot_MW;    % net export
                Pwel_inj_MW(i) = max(Pwel_inj_MW(i), 0);      % no net import in this model (green hydrogen)
                Pwel_inj_MW(i) = min(Pwel_inj_MW(i), Pgrid_con_MW); % grid connection cap
                % Rewrite stochastic flexibility limits
                mpc_w.wind_elec.droop_up(i, 1) = min(Pw_avail_MW,Pgrid_con_MW);
                mpc_w.wind_elec.droop_dwn(i, 1) = max(-Pel_cap_MW+Pel_hot_MW,-Pgrid_con_MW);
            end
        
            % write net injections into the gen table
            mpc_w.gen(wel_idx, PG)   = Pwel_inj_MW;
            mpc_w.gen(wel_idx, PMAX) = Pwel_inj_MW;
            mpc_w.gen(wel_idx, PMIN) = 0;
            
            % j = find(wel_idx == i, 1, 'first');      % WEL unit index 
            

        
            % optional debug prints
            if isfield(opts,'debug') && opts.debug
                fprintf('WEL: mean P_inj=%.2f MW, max=%.2f MW\n', mean(Pwel_inj_MW), max(Pwel_inj_MW));
            end
        end

        % ---- (1b) sample load multiplier + build sampled mpc ----
        [mpc_s, mult, Psys_MW] = sample_loads(mpc_w, ldc, sampler);
        mult_samples(it) = mult;

        % % ---- (1) sample distributed loads + build sampled mpc ----
        % [mpc_s, PD_s_MW, info] = sample_loads_distributed(mpc, ldc, sampler);
        % 
        % % store global draw (from LDC)
        % mult_samples(it) = info.m_global;
        % 
        % % store realized system ratio (should equal m_global if keep_sum=true)
        % mult_realized(it) = info.Psys_MW / info.Psys0_MW;
        % mult=mult_realized(it);
        % 
        % % optional: store the actual system load in MW
        % Psys_samples(it) = info.Psys_MW;
        % Psys_MW = Psys_samples(it);

        % ---- (2) solve CSOPF on sampled case ----
        results = run_csopf_dc(mpc_s, mpopt);

        if ~isfield(results,'success') || results.success == 0
            success(it) = false;
            continue;  % skip failed samples (or count them separately)
        end
        success(it) = true;

        % ---- (3) compute droop-stage congestion index ----
        CI = csopf_droop_CI(results, mpc_s, opts.mode);

        % scalar per-sample metric = contingency-weighted expected CI
        x = CI.E_CI;   % = sum(p_k * CI_system(k))
        E_samples(it) = x;

        % ---- init vector means on first valid sample ----
        if ~have_vec
            mean_RDp = zeros(size(CI.E_RDp));
            mean_RDm = zeros(size(CI.E_RDm));
            mean_LS  = zeros(size(CI.E_LS));
            mean_LFp = zeros(size(CI.E_LFp));
            mean_LFm = zeros(size(CI.E_LFm));
        
            % store indices for interpreting rows
            RD_gen_idx = CI.RD_gen_idx(:);
            LF_bus_idx = CI.LF_lod_idx(:);
        
            have_vec = true;
        end

        % ---- update running mean/variance (Welford) ----
        n = n + 1;
        delta = x - meanE;
        meanE = meanE + delta / n;
        M2 = M2 + delta * (x - meanE);
        % ---- update running means for vector expectations ----
        mean_RDp = mean_RDp + (CI.E_RDp - mean_RDp) / n;
        mean_RDm = mean_RDm + (CI.E_RDm - mean_RDm) / n;
        mean_LS  = mean_LS  + (CI.E_LS  - mean_LS ) / n;
        mean_LFp = mean_LFp + (CI.E_LFp - mean_LFp) / n;
        mean_LFm = mean_LFm + (CI.E_LFm - mean_LFm) / n;

        if n >= 2
            s2 = M2 / (n - 1);
            se = sqrt(s2 / n);
            half = z * se;
        else
            half = Inf;
            se = Inf;
        end

        % ---- optional progress print ----
        if mod(it, opts.print_every) == 0
            fprintf('MC it=%d (valid=%d): meanE=%.4f, halfCI=%.4g, mult=%.3f, Psys=%.1f MW\n', ...
                it, n, meanE, half, mult, Psys_MW);
        end

        % ---- stopping criterion ----
        if n >= opts.Nmin
            rel_ok = (half <= opts.rel_tol * max(abs(meanE), 1e-9));
            abs_ok = (half <= opts.abs_tol);
            if rel_ok || abs_ok
                fprintf('STOP: converged at it=%d (valid=%d): meanE=%.6f, halfCI=%.3g\n', it, n, meanE, half);
                break;
            end
        end
    end

    % pack outputs
    MC = struct();
    MC.meanE = meanE;
    MC.n_valid = n;
    if n >= 2
        MC.varE = M2 / (n - 1);
        MC.se = sqrt(MC.varE / n);
        MC.halfCI = z * MC.se;
    else
        MC.varE = NaN;
        MC.se = NaN;
        MC.halfCI = NaN;
    end
    MC.E_samples = E_samples(1:it);
    MC.mult_samples = mult_samples(1:it);
    MC.success = success(1:it);
    % MC-estimated means of the contingency-expected vectors
    MC.E_RDp = mean_RDp;   % (ngr x 1) in MW
    MC.E_RDm = mean_RDm;   % (ngr x 1) in MW
    MC.E_LS  = mean_LS;    % (nb  x 1) in MW
    MC.E_LFp = mean_LFp;   % (nf  x 1) in MW
    MC.E_LFm = mean_LFm;   % (nf  x 1) in MW
    
    % indices to interpret those vectors
    MC.RD_gen_idx = RD_gen_idx;  % generator rows corresponding to E_RDp/E_RDm
    MC.LF_lod_idx = LF_bus_idx;  % bus indices corresponding to E_LFp/E_LFm
    MC.opts = opts;
end
