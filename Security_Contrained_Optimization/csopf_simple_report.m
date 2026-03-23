function csopf_simple_report(results, mpc)
%CSOPF_SIMPLE_REPORT  Quick text summary of CS-DC-OPF results, including droop stage + CI.
%
%   csopf_simple_report(results, mpc)
%
%   Assumes run_csopf_dc(mpc, mpopt) was used and that add_csopf_dc()
%   created tertiary vars named:
%       'Pg_RDp_k', 'Pg_RDm_k', 'P_LS_k', 'P_LFp_k', 'P_LFm_k'
%   for each contingency k.
%
%   Uses mpc.csopf.Gf, mpc.csopf.R, mpc.csopf.cont(k).dPg_idx, dPd, etc.
%   in *external* indexing.

    define_constants;

    if ~results.success
        fprintf('\n=== CS-DC-OPF FAILED (success = 0) ===\n');
        return;
    end

    baseMVA = results.baseMVA;
    bus     = results.bus;      % external order
    gen     = results.gen;      % external order
    branch0 = mpc.branch;       % base-case branch data (external)
    bus0    = mpc.bus;          % base-case bus data (external)

    nb = size(bus, 1);
    ng = size(gen, 1);
    nl = size(branch0, 1);

    cs   = mpc.csopf;
    cont = cs.cont;
    nc   = numel(cont);

    Gf = cs.Gf(:);    % droop gens (external indices)
    R  = cs.R(:);     % droop coefficients for those gens

    Lf = cs.Lf(:);    % droop loads (may be empty)
    RL = cs.RL(:);

    fprintf('\n============================================\n');
    fprintf('         CS-DC-OPF SIMPLE REPORT\n');
    fprintf('============================================\n');

    %% ---------- Base case summary ----------
    fprintf('\n--- Base case ---\n');
    fprintf('Objective value f = %.4f\n', results.f);

    % base-case Pg (pu & MW)
    Pg0_pu = results.var.val.Pg;          % ng x 1, pu
    Pg0_MW = Pg0_pu * baseMVA;

    Pd_MW = bus(:, PD);                  % MW
    Pgen_tot = sum(Pg0_MW);
    Pload_tot = sum(Pd_MW);

    fprintf('Total generation : %.2f MW\n', Pgen_tot);
    fprintf('Total load       : %.2f MW\n', Pload_tot);
    fprintf('Mismatch (G-L)   : %.4f MW\n', Pgen_tot - Pload_tot);

    % print generator table
    fprintf('\nBase-case generator outputs:\n');
    fprintf('   Gen  Bus   Pg [MW]\n');
    fprintf('  ---------------------\n');
    for g = 1:ng
        fprintf('   %2d   %3d   %8.2f\n', g, gen(g, GEN_BUS), Pg0_MW(g));
    end

    %% ---------- Per-contingency details ----------
    Pd_pu = Pd_MW / baseMVA;    % not strictly needed, but handy

    % base bus-gen incidence (external)
    Bgen = sparse(gen(:, GEN_BUS), (1:ng)', 1, nb, ng);

    for k = 1:nc
        ck = cont(k);

        name_RDp = sprintf('Pg_RDp_%d', k);
        name_RDm = sprintf('Pg_RDm_%d', k);
        name_LS  = sprintf('P_LS_%d',  k);
        name_LFp = sprintf('P_LFp_%d', k);
        name_LFm = sprintf('P_LFm_%d', k);

        fprintf('\n============================================\n');
        fprintf(' Contingency %d\n', k);
        fprintf('   type = %s, prob = %.4f\n', ck.type, ck.prob);
        if strcmpi(ck.type, 'line')
            fprintf('   line outage = %d\n', ck.line);
        elseif strcmpi(ck.type, 'inj')
            if isfield(ck, 'dPg_idx') && ~isempty(ck.dPg_idx)
                fprintf('   gen outage  = %d\n', ck.dPg_idx);
            end
            if isfield(ck, 'dPd') && ~isempty(ck.dPd)
                fprintf('   has dPd profile\n');
            end
        end

        %% ---------- Secondary stage (droop) reconstruction ----------
        % Build ΔP_sys^(c) = a_c' * Pg0_pu + const_c, in pu
        a_c     = zeros(ng, 1);     % external gen indexing
        const_c = 0;
        gtrip   = [];

        if strcmpi(ck.type, 'inj')
            % generator outage
            if isfield(ck, 'dPg_idx') && ~isempty(ck.dPg_idx)
                gtrip = ck.dPg_idx;   % external index
                a_c(gtrip) = 1;       % ΔP_sys depends on P_gtrip
            end
            % load change vector dPd (MW) → pu
            if isfield(ck, 'dPd') && ~isempty(ck.dPd)
                const_c = const_c + sum(ck.dPd) / baseMVA;  % Σ dPd (pu)
            end
        end

        % droop participants: gens + loads (if any)
        Gf_k  = Gf(~ismember(Gf, gtrip));   % droop gens excluding tripped one
        R_k   = R(~ismember(Gf, gtrip));
        ngf_k = numel(Gf_k);

        Lf_k  = Lf;                         % droop load buses (external)
        RL_k  = RL;
        nlf_k = numel(Lf_k);

        % Build combined weights wG, wL and shares alpha (gens) and beta (loads)
        if ngf_k > 0
            wG = 1 ./ R_k;
        else
            wG = 0;
        end
        if nlf_k > 0
            wL = 1 ./ RL_k;
        else
            wL = 0;
        end

        Wtot = sum(wG) + sum(wL);
        if Wtot == 0
            alphaG_k = zeros(ngf_k, 1);
            betaL_k  = zeros(nlf_k, 1);
        else
            alphaG_k = wG / Wtot;       % droop shares for generators
            betaL_k  = wL / Wtot;       % droop shares for loads
        end

        % System imbalance in pu (based on optimal Pg0)
        dP_sys_pu = a_c' * Pg0_pu + const_c;    % scalar pu
        dP_sys_MW = dP_sys_pu * baseMVA;        % same in MW

        %---------- Secondary-stage generator outputs Pg2 (pu/MW) ----------
        Pg2_pu = Pg0_pu;    % start from base case (pu)
        if ~isempty(gtrip)
            Pg2_pu(gtrip) = 0;    % tripped generator output → 0
        end
        for r = 1:ngf_k
            gi = Gf_k(r);             % external gen index
            Pg2_pu(gi) = Pg0_pu(gi) + alphaG_k(r) * dP_sys_pu;
        end
        Pg2_MW = Pg2_pu * baseMVA;

        % print droop table for relevant generators (droop + tripped)
        gens_to_show = unique([Gf_k(:); gtrip(:)]);
        gens_to_show(gens_to_show == 0) = [];
        gens_to_show = sort(gens_to_show);

        Pmin_pu = mpc.gen(:, PMIN) / baseMVA;
        Pmax_pu = mpc.gen(:, PMAX) / baseMVA;

        fprintf('\n--- Secondary stage (droop) generator outputs ---\n');
        if isempty(gens_to_show)
            fprintf('   (No droop generators defined for this contingency)\n');
        else
            fprintf('   Gen   Bus   Pg0_pu   Pg2_pu   Pmin_pu   Pmax_pu   Above_Pmin  Below_Pmax\n');
            fprintf('  -------------------------------------------------------------------------\n');
            for gg = gens_to_show'
                above_min  = Pg2_pu(gg) >= Pmin_pu(gg) - 1e-6;
                below_max  = Pg2_pu(gg) <= Pmax_pu(gg) + 1e-6;
                fprintf('   %3d  %4d  %7.4f  %7.4f  %8.4f  %8.4f     %5s       %5s\n', ...
                    gg, gen(gg, GEN_BUS), Pg0_pu(gg), Pg2_pu(gg), ...
                    Pmin_pu(gg), Pmax_pu(gg), ...
                    tf_str(above_min), tf_str(below_max));
            end
        end

        %---------- Secondary-stage frequency-sensitive loads PD2 (pu/MW) ----------
        fprintf('\n--- Secondary stage (droop) frequency-sensitive loads ---\n');
        if nlf_k == 0
            fprintf('   (No frequency-sensitive loads defined)\n');
        else
            fprintf('   Bus   PD0_pu   PD2_pu\n');
            fprintf('  ------------------------\n');

            % Start from base loads (MW) + explicit dPd (if any)
            Pd_c_MW = Pd_MW;
            if strcmpi(ck.type, 'inj') && isfield(ck, 'dPd') && ~isempty(ck.dPd)
                Pd_c_MW = Pd_c_MW + ck.dPd(:);    % MW
            end

            for r = 1:nlf_k
                bj       = Lf_k(r);                 % bus index (external)
                PD0_MW   = Pd_MW(bj);               % base MW
                PD0_pu   = PD0_MW / baseMVA;

                % droop response: PD2 = PD0 + dPd_j - beta * ΔP_sys
                PD2_MW   = Pd_c_MW(bj) - betaL_k(r) * dP_sys_MW;
                PD2_pu   = PD2_MW / baseMVA;

                fprintf('   %3d  %7.4f  %7.4f\n', bj, PD0_pu, PD2_pu);
            end
        end

        %% ---------- Droop-stage congestion index CI (with load droop) ----------
        % Build contingency network for PTDF
        branch_c = branch0;
        if strcmpi(ck.type, 'line')
            branch_c(ck.line, BR_STATUS) = 0;
        end
        PTDF_c = makePTDF(baseMVA, bus0, branch_c);   % nl x nb

        % Rebuild Pd_c_MW including dPd and load droop (reusing betaL_k)
        Pd_c_MW = Pd_MW;
        if strcmpi(ck.type, 'inj') && isfield(ck, 'dPd') && ~isempty(ck.dPd)
            Pd_c_MW = Pd_c_MW + ck.dPd(:);           % explicit load change
        end
        if nlf_k > 0 && any(betaL_k)
            for r = 1:nlf_k
                bj = Lf_k(r);
                Pd_c_MW(bj) = Pd_c_MW(bj) - betaL_k(r) * dP_sys_MW;
            end
        end

        % Now droop-stage bus injections (MW)
        % Generators at Pg2_MW, loads at Pd_c_MW (with dPd + droop)
        Bgen_c = Bgen;
        if ~isempty(gtrip)
            Bgen_c(:, gtrip) = 0;
        end
        Pinj_droop_MW = Bgen_c * Pg2_MW - Pd_c_MW;

        % Line flows under droop-only response
        F_droop_MW = PTDF_c * Pinj_droop_MW;
        % save('F_droop_0.001_FSM.mat', 'F_droop_MW');

        % Line ratings (MW) for contingency network
        rateC_c = branch_c(:, RATE_C);
        rateC_c(rateC_c == 0) = Inf;

        CI_droop_line   = abs(F_droop_MW) ./ rateC_c;
        CI_droop_system = mean(CI_droop_line);
        [CImax_droop, lidx_droop] = max(CI_droop_line);

        fprintf('\n--- Droop-stage congestion index ---\n');
        fprintf('   mean_CI_system = %.4f\n', CI_droop_system);
        if CImax_droop > 0 && ~isempty(lidx_droop)
            fprintf('   Worst line (droop) %d: |F| = %.2f MW, Fmax = %.2f MW, CI = %.4f\n', ...
                lidx_droop, F_droop_MW(lidx_droop), rateC_c(lidx_droop), CImax_droop);
        end

        %% ---------- Tertiary stage summary (if vars exist) ----------
        if ~isfield(results.var.val, name_RDp)
            fprintf('\n   (No tertiary vars found for this contingency)\n');
            continue;
        end

        RDp_pu = results.var.val.(name_RDp);   % ng x 1
        RDm_pu = results.var.val.(name_RDm);   % ng x 1
        LS_pu  = results.var.val.(name_LS);    % nb x 1
        LFp_pu = results.var.val.(name_LFp);   % nf x 1
        LFm_pu = results.var.val.(name_LFm);   % nf x 1

        RDp_MW = RDp_pu * baseMVA;
        RDm_MW = RDm_pu * baseMVA;
        LS_MW  = LS_pu  * baseMVA;
        LFp_MW = LFp_pu * baseMVA;
        LFm_MW = LFm_pu * baseMVA;

        % Sums
        sum_RDp = sum(RDp_MW);
        sum_RDm = sum(RDm_MW);
        sum_LS  = sum(LS_MW);
        sum_LFp = sum(LFp_MW);
        sum_LFm = sum(LFm_MW);

        fprintf('\n--- Tertiary actions (MW) ---\n');
        fprintf('   sum(RDp) = %10.4f MW\n', sum_RDp);
        fprintf('   sum(RDm) = %10.4f MW\n', sum_RDm);
        fprintf('   sum(LS)  = %10.4f MW\n', sum_LS);
        fprintf('   sum(LFp) = %10.4f MW\n', sum_LFp);
        fprintf('   sum(LFm) = %10.4f MW\n', sum_LFm);

        % Net tertiary adjustment (should be ~0 if balance constraint is right):
        net_ter = sum_RDp - sum_RDm - sum_LS - sum_LFp + sum_LFm;
        fprintf('   Net tertiary adj (RDp - RDm - LS - LFp + LFm) = %10.4f MW\n', net_ter);

        %% non-zero LS per bus
        tol = 1e-3;  % MW
        idx_LS = find(abs(LS_MW) > tol);
        if ~isempty(idx_LS)
            fprintf('\n   Non-zero load shedding:\n');
            fprintf('      Bus    LS [MW]\n');
            for ii = 1:length(idx_LS)
                b = idx_LS(ii);
                fprintf('      %3d   %8.3f\n', bus(b, BUS_I), LS_MW(b));
            end
        else
            fprintf('\n   No load shedding (LS) in this contingency.\n');
        end

        %% non-zero flex load changes (LFp, LFm)
        flex_buses = mpc.csopf.flex_buses(:);
        nf = numel(flex_buses);

        if nf > 0
            LFp_bus_MW = zeros(nb, 1);
            LFm_bus_MW = zeros(nb, 1);
            LFp_bus_MW(flex_buses) = LFp_MW;
            LFm_bus_MW(flex_buses) = LFm_MW;

            idx_LF = find(abs(LFp_bus_MW - LFm_bus_MW) > tol);
            if ~isempty(idx_LF)
                fprintf('\n   Non-zero flexible load changes (LFp/LFm):\n');
                fprintf('      Bus    LFp [MW]    LFm [MW]\n');
                for ii = 1:length(idx_LF)
                    b = idx_LF(ii);
                    fprintf('      %3d   %8.3f   %8.3f\n', ...
                        bus(b, BUS_I), LFp_bus_MW(b), LFm_bus_MW(b));
                end
            else
                fprintf('\n   No flexible load changes (LF) in this contingency.\n');
            end
        end
    end

    fprintf('\n============================================\n');
    fprintf(' End of CS-DC-OPF simple report\n');
    fprintf('============================================\n\n');
end

%% helper: boolean → 'true'/'false'
function s = tf_str(tf)
    if tf
        s = 'true';
    else
        s = 'false';
    end
end
