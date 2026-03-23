function CI = csopf_droop_CI(results, mpc, mode)
%CSOPF_DROOP_CI  Compute droop-stage congestion indices for CSOPF results.
%
%   CI = csopf_droop_CI(results, mpc, mode)
%
% Inputs
%   results : MATPOWER results struct from run_csopf_dc (external indexing)
%   mpc     : MATPOWER case used for the solve (same loads as results)
%   mode    : (optional) 'mean' (default) or 'max'
%             - 'mean' => CI_system = mean(|F|/Fmax) across lines
%             - 'max'  => CI_system = max(|F|/Fmax) across lines (classic)
%
% Output struct CI
%   .CI_system(k)      : system congestion index per contingency
%   .CI_line(:,k)      : line congestion indices (nl x nc)
%   .worst_line(k)     : line index giving maximum CI
%   .CImax(k)          : maximum CI value per contingency
%   .F_droop_MW(:,k)   : droop-stage line flows (MW), nl x nc
%   .E_CI              : expected system CI = sum(p_k * CI_system(k))
%   .prob(:)           : contingency probabilities
%
% Notes
% - Uses RATE_C for limits (as you prefer).
% - Includes optional frequency-sensitive load droop if mpc.csopf.Lf non-empty.
% - Assumes DC model and Pg values are in MW in results.gen(:,PG).

    define_constants;

    if nargin < 3 || isempty(mode)
        mode = 'mean';
    end

    if ~isfield(mpc, 'csopf') || ~isfield(mpc.csopf, 'cont')
        error('csopf_droop_CI: mpc.csopf.cont not found.');
    end
    cont = mpc.csopf.cont;
    nc   = numel(cont);

    baseMVA = mpc.baseMVA;
    bus0    = mpc.bus;
    branch0 = mpc.branch;
    gen0    = mpc.gen;

    nb = size(bus0,1);
    nl = size(branch0,1);
    ng = size(gen0,1);

    % Base loads (MW) in external bus order
    Pd_MW = bus0(:, PD);

    % Generator-bus incidence (nb x ng) in external indexing
    Bgen = sparse(gen0(:, GEN_BUS), (1:ng)', 1, nb, ng);

    % Precompute PTDF for base topology (reused for non-line contingencies)
    PTDF_base = makePTDF(baseMVA, bus0, branch0);

    % Droop configuration
    cs = mpc.csopf;
    Gf = cs.Gf(:);
    R  = cs.R(:);

    Lf = [];
    RL = [];
    if isfield(cs, 'Lf') && ~isempty(cs.Lf)
        Lf = cs.Lf(:);
        RL = cs.RL(:);
    end

    % Base dispatch Pg0 from results (MW, external gen order)
    Pg0_MW = results.gen(:, PG);
    Pg0_pu = Pg0_MW / baseMVA;

    % Redispatch outputs
    redispatch_gen_idx = mpc.csopf.redispatch_gen_idx(:);
    flex_buses = mpc.csopf.flex_buses(:);
    ngr = numel(redispatch_gen_idx); nf = numel(flex_buses);

    %% Pull Tertiary outputs
    % storage per contingency (PU)
    RDp_pu = zeros(ngr, nc);
    RDm_pu = zeros(ngr, nc);
    LS_pu  = zeros(nb,  nc);
    LFp_pu = zeros(nf,  nc);
    LFm_pu = zeros(nf,  nc);
    % pull variables
    for k = 1:nc
        name_RDp = sprintf('Pg_RDp_%d', k);
        name_RDm = sprintf('Pg_RDm_%d', k);
        name_LS  = sprintf('P_LS_%d',  k);
        name_LFp = sprintf('P_LFp_%d', k);
        name_LFm = sprintf('P_LFm_%d', k);

        if isfield(results.var.val, name_RDp)
            RDp_pu(:,k) = results.var.val.(name_RDp);
        end
        if isfield(results.var.val, name_RDm)
            RDm_pu(:,k) = results.var.val.(name_RDm);
        end

        % LS is nb×1, indexed by bus directly
        if isfield(results.var.val, name_LS)
            LS_pu(:,k) = results.var.val.(name_LS);
        end

        % LFp/LFm are nf×1 aligned with flex_buses list
        if nf > 0
            if isfield(results.var.val, name_LFp)
                LFp_pu(:,k) = results.var.val.(name_LFp);
            end
            if isfield(results.var.val, name_LFm)
                LFm_pu(:,k) = results.var.val.(name_LFm);
            end
        end
    end

    % convert to MW
    RDp_MW = RDp_pu * baseMVA;
    RDm_MW = RDm_pu * baseMVA;
    LS_MW  = LS_pu  * baseMVA;
    LFp_MW = LFp_pu * baseMVA;
    LFm_MW = LFm_pu * baseMVA;

    % Allocate outputs
    CI_line    = zeros(nl, nc);
    CI_system  = zeros(nc, 1);
    worst_line = zeros(nc, 1);
    CImax      = zeros(nc, 1);
    F_droop_MW = zeros(nl, nc);
    prob       = zeros(nc, 1);

    for k = 1:nc
        ck = cont(k);
        prob(k) = ck.prob;

        %% ---- Secondary stage droop reconstruction: Pg2 ----
        a_c     = zeros(ng, 1);   % vector selecting tripped gen in ΔP_sys
        const_c = 0;              % Σ dPd in pu
        gtrip   = [];

        if strcmpi(ck.type, 'inj')
            if isfield(ck, 'dPg_idx') && ~isempty(ck.dPg_idx)
                gtrip = ck.dPg_idx;     % external gen index
                a_c(gtrip) = 1;         % ΔP_sys includes +Pgtrip
            end
            if isfield(ck, 'dPd') && ~isempty(ck.dPd)
                const_c = const_c + sum(ck.dPd(:)) / baseMVA;  % pu
            end
        end

        % droop gens excluding tripped gen
        Gf_k = Gf(~ismember(Gf, gtrip));
        R_k  = R(~ismember(Gf, gtrip));
        ngf_k = numel(Gf_k);

        % FSM loads
        Lf_k  = Lf;
        RL_k  = RL;
        nlf_k = numel(Lf_k);

        % weights
        if ngf_k > 0
            wG = 1 ./ R_k;
        else
            wG = zeros(0,1);
        end
        if nlf_k > 0
            wL = 1 ./ RL_k;
        else
            wL = zeros(0,1);
        end
        Wtot = sum(wG) + sum(wL);

        if Wtot > 0 && ngf_k > 0
            alphaG_k = wG / Wtot;
        else
            alphaG_k = zeros(ngf_k,1);
        end
        if Wtot > 0 && nlf_k > 0
            betaL_k = wL / Wtot;
        else
            betaL_k = zeros(nlf_k,1);
        end

        % ΔP_sys in pu and MW
        dP_sys_pu = a_c' * Pg0_pu + const_c;
        dP_sys_MW = dP_sys_pu * baseMVA;

        % Pg2 (pu)
        Pg2_pu = Pg0_pu;
        if ~isempty(gtrip)
            Pg2_pu(gtrip) = 0;
        end
        for r = 1:ngf_k
            gi = Gf_k(r);
            Pg2_pu(gi) = Pg0_pu(gi) + alphaG_k(r) * dP_sys_pu;
        end
        Pg2_MW = Pg2_pu * baseMVA;

        %% ---- Build contingency PTDF (line outages only change topology) ----
        branch_c = branch0;
        if strcmpi(ck.type, 'line')
            branch_c(ck.line, BR_STATUS) = 0;
            PTDF_c = makePTDF(baseMVA, bus0, branch_c);
        else
            PTDF_c = PTDF_base;
        end

        %% ---- Build Pd_c_MW including explicit dPd and FSM droop ----
        Pd_c_MW = Pd_MW;

        if strcmpi(ck.type, 'inj') && isfield(ck, 'dPd') && ~isempty(ck.dPd)
            Pd_c_MW = Pd_c_MW + ck.dPd(:);
        end

        % Apply FSM droop: PD_j^c = PD_j + dPd_j - beta_j * ΔP_sys
        if nlf_k > 0 && any(betaL_k)
            for r = 1:nlf_k
                bj = Lf_k(r);                    % external bus index
                Pd_c_MW(bj) = Pd_c_MW(bj) - betaL_k(r) * dP_sys_MW;
            end
        end

        %% ---- Bus injections and droop flows ----
        Bgen_c = Bgen;
        if ~isempty(gtrip)
            Bgen_c(:, gtrip) = 0;
        end

        Pinj_droop_MW = Bgen_c * Pg2_MW - Pd_c_MW;    % nb x 1
        Fk_MW = PTDF_c * Pinj_droop_MW;               % nl x 1
        F_droop_MW(:,k) = Fk_MW;

        %% ---- CI ----
        rateC = branch_c(:, RATE_C);
        rateC(rateC == 0) = Inf;

        CI_k = abs(Fk_MW) ./ rateC;
        CI_line(:,k) = CI_k;

        [CImax(k), worst_line(k)] = max(CI_k);

        switch lower(mode)
            case 'mean'
                CI_system(k) = mean(CI_k(isfinite(CI_k))); % ignore Inf-limits
            case 'max'
                CI_system(k) = CImax(k);
            otherwise
                error('csopf_droop_CI: mode must be ''mean'' or ''max''.');
        end
    end

    % Expected value across contingencies
    E_CI = sum(prob .* CI_system);
    E_RDp_gen_MW = RDp_MW * prob;     % ngr×1
    E_RDm_gen_MW = RDm_MW * prob;     % ngr×1
    E_LS_bus_MW  = LS_MW  * prob;     % nb×1
    E_LFp_bus_MW = LFp_MW * prob; % nf×1 
    E_LFm_bus_MW = LFm_MW * prob; % nf×1

    % Pack output
    CI = struct();
    CI.CI_system  = CI_system;
    CI.E_RDp=E_RDp_gen_MW;
    CI.E_RDm=E_RDm_gen_MW;
    CI.RD_gen_idx=redispatch_gen_idx;

    CI.E_LFp=E_LFp_bus_MW;
    CI.E_LFm=E_LFm_bus_MW;
    CI.LF_lod_idx=flex_buses;

    CI.E_LS=E_LS_bus_MW;
    CI.CI_line    = CI_line;
    CI.worst_line = worst_line;
    CI.CImax      = CImax;
    CI.F_droop_MW  = F_droop_MW;
    CI.prob       = prob;
    CI.E_CI       = E_CI;
    CI.mode       = mode;
end
