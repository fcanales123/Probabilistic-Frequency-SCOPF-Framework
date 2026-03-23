function om = add_csopf_dc(om, mpc)
    define_constants;
    baseMVA = mpc.baseMVA;
    bus     = mpc.bus;
    branch  = mpc.branch;
    gen     = mpc.gen;

    nb = size(bus, 1);
    nl = size(branch, 1);
    ng = size(gen, 1);

    %% ----- CS-OPF configuration -----
    cs           = mpc.csopf;
    wel          = mpc.wind_elec;

    % If ext2int created an 'order' field, use it
    if isfield(mpc, 'order')
        o = mpc.order;

        % ---- generator-related indices ----
        if ~isempty(cs.Gf)
            % external gen indices -> internal gen row indices
            cs.Gf = o.gen.e2i(cs.Gf);
        end
        for kk = 1:numel(cs.cont) % Generator outages
            if isfield(cs.cont(kk), 'dPg_idx') && ~isempty(cs.cont(kk).dPg_idx)
                cs.cont(kk).dPg_idx = o.gen.e2i(cs.cont(kk).dPg_idx);
            end
        end

        % ---- bus-related indices ----
        if ~isempty(cs.flex_buses)
            cs.flex_buses = o.bus.e2i(cs.flex_buses);
        end
        if ~isempty(cs.Lf)
            cs.Lf = o.bus.e2i(cs.Lf);
        end

        % dPd is a bus-indexed vector; reorder it to internal bus order
        for kk = 1:numel(cs.cont)
            if isfield(cs.cont(kk), 'dPd') && ~isempty(cs.cont(kk).dPd)
                dPd_ext = cs.cont(kk).dPd(:);   % external bus order
                dPd_int = zeros(nb,1);
                for b_ext = 1:length(o.bus.e2i)
                    b_int = o.bus.e2i(b_ext);
                    if b_int > 0
                        dPd_int(b_int) = dPd_ext(b_ext);
                    end
                end
                cs.cont(kk).dPd = dPd_int;  % overwrite with internal ordering
            end
        end
    end

    Gf           = cs.Gf(:);          % droop gens indices
    wel_id       = wel.idx;           % wel droop gens indices
    R            = cs.R(:);           % droop coefficients for those gens
    C_RD         = cs.C_RD(:);        % redispatch cost (ng x 1)
    C_LS         = cs.C_LS(:);        % load shedding cost (nb x 1)
    C_LF         = cs.C_LF(:);        % load flexibility cost (nb x 1)
    flex_buses   = cs.flex_buses(:);  % buses with flexible loads
    redispatch_gens   = cs.redispatch_gen_idx(:);  % redispatchable gens
    PLF_max      = cs.PLF_max(:);     % |ΔP_LF|max per flex bus
    cont         = cs.cont;           % contingencies struct array
    nc           = numel(cont);

    Lf      = cs.Lf(:);     % bus indices of droop loads
    RL      = cs.RL(:);     % their droop R
    PD_min  = cs.PD_min(:); % MW, aligned with Lf
    PD_max  = cs.PD_max(:); % MW, aligned with Lf

    %% ----- PTDF for base network (reuse for non-topology contingencies) -----
    PTDF_base = makePTDF(baseMVA, bus, branch);

    %% ----- Get Pg variable indices in the opf_model -----
    idx = om.get_idx('var');      % get all index info
    
    Pg_i1 = idx.i1.Pg;   % starting index of Pg block in x
    Pg_iN = idx.iN.Pg;   % ending index of Pg block in x
    Pg_N  = Pg_iN - Pg_i1 + 1;
    
    assert(Pg_N == ng, 'Pg var size mismatch with gen matrix');
    assert(Pg_N == size(mpc.gen,1), 'Mismatch: Pg vars != number of generators');

    %% ----- Bus-generator incidence matrix B (nb x ng) -----
    B = sparse(gen(:, GEN_BUS), (1:ng)', 1, nb, ng);

    %% ----- Precompute base load vector (pu) -----
    Pd = bus(:, PD) / baseMVA;

    %% ----- Indices for redispatchable generators -----
    nrd = numel(redispatch_gens);

    %% ----- Flexible-load incidence matrix Sflex (nb x nf) -----
    nf     = numel(flex_buses);
    Sflex  = sparse(flex_buses, (1:nf)', 1, nb, nf);

    %% ----- Global droop participation factors (only used for stage-2 gen limits) -----
    invR   = 1 ./ R;
    alpha  = invR / sum(invR);   % |Gf| x 1 (global, will be filtered per cont.)

    %% ----- Loop over contingencies -----
    for k = 1:nc
        ck = cont(k);
        pc = ck.prob;    % probability of contingency k (for weighting costs)

        %% =========================
        %  (1) Add tertiary variables
        %% =========================

        name_RDp = sprintf('Pg_RDp_%d', k);
        name_RDm = sprintf('Pg_RDm_%d', k);
        name_LS  = sprintf('P_LS_%d',  k);
        name_LFp = sprintf('P_LFp_%d', k);
        name_LFm = sprintf('P_LFm_%d', k);

        % Redispatch variables only for nrd generators
        om.add_var(name_RDp, nrd, 0, zeros(nrd,1), Inf(nrd,1));
        om.add_var(name_RDm, nrd, 0, zeros(nrd,1), Inf(nrd,1));

        % Load shedding (nb)
        om.add_var(name_LS, nb, 0, zeros(nb,1), Inf(nb,1));

        % Flexible load +/- (nf each)
        om.add_var(name_LFp, nf, 0, zeros(nf,1), PLF_max);
        om.add_var(name_LFm, nf, 0, zeros(nf,1), PLF_max);

        %% =========================
        %  (2) Add tertiary costs
        %% =========================

        % Redispatch + (Pg_RDp_k)
        Q_RD = sparse(nrd, nrd);                 % no quadratic term
        c_RD = pc * C_RD(:);                   % ng x 1
        om.add_quad_cost( ...
            sprintf('cost_RDp_%d', k), ...     % unique name per contingency
            Q_RD, c_RD, 0, {name_RDp} );

        % Redispatch - (Pg_RDm_k)
        om.add_quad_cost( ...
            sprintf('cost_RDm_%d', k), ...
            Q_RD, c_RD, 0, {name_RDm});

        % Load shedding (P_LS_k), all buses
        Q_LS = sparse(nb, nb);
        c_LS = pc * C_LS(:);                   % nb x 1
        om.add_quad_cost( ...
            sprintf('cost_LS_%d', k), ...
            Q_LS, c_LS, 0, {name_LS} );

        % Flexible load +/- (only on flex_buses)
        Q_LF     = sparse(nf, nf);
        % CLF_flex = C_LF(flex_buses);           % nf x 1
        c_LF     = pc * C_LF(:);

        om.add_quad_cost( ...
            sprintf('cost_LFp_%d', k), ...
            Q_LF, c_LF, 0, {name_LFp} );

        om.add_quad_cost( ...
            sprintf('cost_LFm_%d', k), ...
            Q_LF, c_LF, 0, {name_LFm} );

        %% =========================
        %  (3) Stage-2 droop-based generator limits (no RD here)
        %% =========================

        %----- compute ΔP_sys^(c) = a_c' * Pg + const_c (symbolic) -----
        a_c     = zeros(ng,1);
        const_c = 0;
        gtrip   = [];

        if strcmpi(ck.type, 'inj')
            % Example: generator outage
            if isfield(ck, 'dPg_idx') && ~isempty(ck.dPg_idx)
                gtrip = ck.dPg_idx;
                a_c(gtrip) = 1;   % ΔP_sys = P_gtrip
            end
            if isfield(ck, 'dPd') && ~isempty(ck.dPd)
                const_c = const_c + sum(ck.dPd)/baseMVA;  % Σ ΔPd (pu)
            end
        end

        %===== DROOP PARTICIPATION (gens + loads) =====
        % Gens for this contingency
        Gf_k  = Gf(~ismember(Gf, gtrip));
        R_k   = R(~ismember(Gf, gtrip));
        ngf_k = numel(Gf_k);

        % Loads (we assume all in Lf participate in all contingencies;
        % you can filter per contingency if needed)
        Lf_k  = Lf;
        RL_k  = RL;
        nlf_k = numel(Lf_k);

        % Build combined weights
        wG = 1 ./ R_k;    % |Gf_k| x 1
        wL = 1 ./ RL_k;   % |Lf_k| x 1
        Wtot = sum(wG) + sum(wL);

        alphaG_k = wG / Wtot;   % droop shares for gens
        betaL_k  = wL / Wtot;   % droop shares for loads

        %----- generator droop limits (same as before, just use alphaG_k) -----
        A_Pg_gen  = sparse(ngf_k, Pg_N);
        l_gen = zeros(ngf_k, 1);
        u_gen = zeros(ngf_k, 1);

        for r = 1:ngf_k
            gi      = Gf_k(r);          % gen index in gen matrix
            alpha_r = alphaG_k(r);      % droop share

            if ismember(gi, wel_id)
                Pmax_i_pu = wel.droop_up(wel_id==gi) / baseMVA;
                Pmin_i_pu = wel.droop_dwn(wel_id==gi) / baseMVA;
            else
                Pmax_i_pu = gen(gi, PMAX) / baseMVA;
                Pmin_i_pu = gen(gi, PMIN) / baseMVA;
            end

            % row = e_i' + alpha_r * a_c'
            A_Pg_gen(r, :)  = alpha_r * a_c.';
            A_Pg_gen(r, gi) = A_Pg_gen(r, gi) + 1;

            l_gen(r) = Pmin_i_pu - alpha_r * const_c;
            u_gen(r) = Pmax_i_pu - alpha_r * const_c;
        end

        if ngf_k > 0
            om.add_lin_constraint( ...
                sprintf('Pg_droop_cont_%d', k), ...
                A_Pg_gen, l_gen, u_gen, ...
                {'Pg'} );
        end

        %===== load droop limits (new) =====
        % For each droop load at bus index bj = Lf_k(r),
        % enforce PD_min <= PD_j^c <= PD_max via Pg (using ΔP_sys)
        %
        % PD_j^c = PD_j - beta_r * (a_c' * Pg + const_c)
        % all quantities in pu

        nlf_k = numel(Lf_k);
        A_Pg_load = sparse(2*nlf_k, Pg_N);   % 2 rows per droop load
        l_load    = -Inf(2*nlf_k, 1);        % one-sided upper bounds
        u_load    = zeros(2*nlf_k, 1);

        for r = 1:nlf_k
            bj     = Lf_k(r);       % bus index (internal ordering)
            beta_r = betaL_k(r);   % droop share for this load

            % Base demand and bounds in pu
            PD_j_MW   = bus(bj, PD);
            PD_min_j  = PD_min(r);     % MW, aligned with Lf
            PD_max_j  = PD_max(r);     % MW, aligned with Lf

            PD_j_pu     = PD_j_MW  / baseMVA;
            PD_min_j_pu = PD_min_j / baseMVA;
            PD_max_j_pu = PD_max_j / baseMVA;

            % Row indices in A_Pg_load
            rowL = 2*r - 1;   % lower bound on PD^c
            rowU = 2*r;       % upper bound on PD^c

            % --- Lower bound: PD_min <= PD_j - beta_r * (a_c'Pg + const_c)
            % => beta_r * a_c' * Pg <= PD_j_pu - PD_min_j_pu - beta_r * const_c
            A_Pg_load(rowL, :) =  beta_r * a_c.';   % 1 x ng
            u_load(rowL)       =  PD_j_pu - PD_min_j_pu - beta_r * const_c;

            % --- Upper bound: PD_j - beta_r * (a_c'Pg + const_c) <= PD_max
            % => -beta_r * a_c' * Pg <= PD_max_j_pu - PD_j_pu + beta_r * const_c
            A_Pg_load(rowU, :) = -beta_r * a_c.';   % 1 x ng
            u_load(rowU)       =  PD_max_j_pu - PD_j_pu + beta_r * const_c;
        end

        if nlf_k > 0
            om.add_lin_constraint( ...
                sprintf('PD_droop_cont_%d', k), ...
                A_Pg_load, l_load, u_load, ...
                {'Pg'} );
        end

        %% =========================
        %  (4) Tertiary-stage network / flow constraints (all in p.u.)
        %% =========================
        
        % Build contingency network for PTDF
        bus_c    = bus;
        branch_c = branch;
        
        if strcmpi(ck.type, 'line')
            % line outage
            branch_c(ck.line, BR_STATUS) = 0;
            PTDF_c = makePTDF(baseMVA, bus_c, branch_c);
        else
            % injection contingency: same topology
            PTDF_c = PTDF_base;
        end
        
        % Bus-generator incidence for contingency (zero out tripped gen)
        B_c = B;
        if ~isempty(gtrip)
            B_c(:, gtrip) = 0;
        end
        
        % Tertiary-stage injection mapping in p.u.:
        % P_inj^c(pu) = B_c * Pg
        %             + B_c * RDp_k - B_c * RDm_k
        %             - (Pd + dPd_pu)
        %             + LS_k
        %             - Sflex*LFp_k + Sflex*LFm_k
        %
        % We encode this as:
        %   P_inj^c = A_total * x + const_inj
        % where x = [Pg; RDp_k; RDm_k; LS_k; LFp_k; LFm_k]
        
        A_Pg_inj  = B_c;         % nb x ng, Pg in pu
        A_RD_inj_k = B_c(:, redispatch_gens);   % nb x nrd (THE KEY CHANGE)

        A_RDp_inj = A_RD_inj_k;         % nb x ng, RDp in pu
        A_RDm_inj = -A_RD_inj_k;        % nb x ng, RDm in pu
        A_LS_inj  = speye(nb);   % nb x nb, LS in pu
        A_LFp_inj = -Sflex;      % nb x nf, LFp in pu
        A_LFm_inj =  Sflex;      % nb x nf, LFm in pu
        
        % Constant injection term in p.u.: base loads + optional dPd
        const_inj = -Pd;         % -(P^D) in pu
        if strcmpi(ck.type, 'inj') && isfield(ck, 'dPd') && ~isempty(ck.dPd)
            const_inj = const_inj - ck.dPd / baseMVA;   % add extra load change in pu
        end

        % ---- System power-balance equality for contingency k ----
        one = ones(1, nb);

        A_balance = [one*A_Pg_inj, one*A_RDp_inj, one*A_RDm_inj, ...
                     one*A_LS_inj, one*A_LFp_inj, one*A_LFm_inj];
        rhs_bal   = - one * const_inj;      % scalar (pu)
        
        om.add_lin_constraint( ...
            sprintf('bal_ter_cont_%d', k), ...
            A_balance, rhs_bal, rhs_bal, ...
            {'Pg', name_RDp, name_RDm, name_LS, name_LFp, name_LFm} );
        
        % Map to flows (p.u.): F_c = PTDF_c * P_inj^c
        A_Pg_flow  = PTDF_c * A_Pg_inj;   % nl x ng
        A_RDp_flow = PTDF_c * A_RDp_inj;  % nl x nrd
        A_RDm_flow = PTDF_c * A_RDm_inj;  % nl x nrd
        A_LS_flow  = PTDF_c * A_LS_inj;   % nl x nb
        A_LFp_flow = PTDF_c * A_LFp_inj;  % nl x nf
        A_LFm_flow = PTDF_c * A_LFm_inj;  % nl x nf
        
        F_const = PTDF_c * const_inj;     % nl x 1, p.u.
        
        % Assemble A_flow over [Pg, RDp_k, RDm_k, LS_k, LFp_k, LFm_k]
        A_flow = [A_Pg_flow, A_RDp_flow, A_RDm_flow, ...
                  A_LS_flow, A_LFp_flow, A_LFm_flow];
        
        % Line limits (p.u.) for contingency network
        rateA = branch_c(:, RATE_A) / baseMVA;    % p.u.
        rateA(rateA == 0) = Inf;                  % 0 => no limit
        
        % We want: -rateA <= F_const + A_flow * x <= rateA
        % → l_flow <= A_flow * x <= u_flow :
        l_flow = -rateA - F_const;
        u_flow =  rateA - F_const;
        
        % Add tertiary flow constraints for contingency k
        om.add_lin_constraint( ...
            sprintf('flow_ter_cont_%d', k), ...
            A_flow, l_flow, u_flow, ...
            {'Pg', name_RDp, name_RDm, name_LS, name_LFp, name_LFm} );

    end
end
