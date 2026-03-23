function [ranked_indices, idx, risk, impact, prob] = rank_contingencies(mpc, top_n)
    define_constants;

    % 1) Base-case DC-OPF
    res_base = rundcopf(mpc, mpoption('verbose', 0));
    if ~res_base.success, error('Base case DC-OPF failed'); end

    nb = size(mpc.bus, 1);
    nl = size(mpc.branch, 1);
    baseMVA = mpc.baseMVA;

    % Base flows/limits in p.u.
    F_base = res_base.branch(:, PF) / baseMVA;
    Limits = mpc.branch(:, RATE_A) / baseMVA;
    Limits(Limits == 0) = 99;

    % Sensitivity factors
    PTDF = makePTDF(baseMVA, mpc.bus, mpc.branch);

    % LODF: use base-case in-service branch statuses
    % (makeLODF expects DC model; PTDF already computed)
    LODF = makeLODF(mpc.branch, PTDF);

    nc = numel(mpc.csopf.cont);
    impact = zeros(nc, 1);
    prob   = zeros(nc, 1);
    risk   = zeros(nc, 1);

    fprintf('Ranking %d contingencies (risk = prob x impact)...\n', nc);

    for k = 1:nc
        ck = mpc.csopf.cont(k);

        % Probability (default 1 if missing)
        if isfield(ck, 'prob') && ~isempty(ck.prob)
            prob(k) = ck.prob;
        else
            prob(k) = 1.0;
        end

        dP_lost = 0;
        tertiary_stress = 0;

        if strcmpi(ck.type, 'inj')
            % --- Generator trip as lost injection ---
            gen_idx = ck.dPg_idx;
            dP_lost = res_base.gen(gen_idx, PG) / baseMVA;

            gen_bus = mpc.gen(gen_idx, GEN_BUS);
            dInj = zeros(nb, 1);
            dInj(gen_bus) = -dP_lost;

            % Post-contingency flow estimate
            F_post = F_base + PTDF * dInj;

            violations = max(0, abs(F_post) - Limits);
            tertiary_stress = sum(violations);

        elseif strcmpi(ck.type, 'line')
            % --- Line outage using LODF ---
            if ~isfield(ck, 'line') || ck.line < 1 || ck.line > nl
                error('Contingency %d has invalid line index.', k);
            end
            ell = ck.line;

            % If the outaged line is out-of-service in base case, skip/zero impact
            if mpc.branch(ell, BR_STATUS) == 0
                impact(k) = 0;
                risk(k)   = 0;
                continue;
            end

            % LODF-based redistribution
            F_post = F_base + LODF(:, ell) * F_base(ell);

            % Ensure outaged line is zero (numerical safety)
            F_post(ell) = 0;

            violations = max(0, abs(F_post) - Limits);
            tertiary_stress = sum(violations);

        else
            % Unknown type: leave at zero (or set a penalty if you prefer)
            dP_lost = 0;
            tertiary_stress = 0;
        end

        % Impact model (your weights)
        impact(k) = dP_lost + (5.0 * tertiary_stress);

        % Risk = probability x impact
        risk(k) = prob(k) * impact(k);
    end

    % Sort by risk
    [~, idx] = sort(risk, 'descend');
    ranked_indices = idx(1:min(top_n, nc));

    fprintf('Top %d most risky contingencies identified.\n', numel(ranked_indices));
end
