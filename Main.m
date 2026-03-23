%% --- Global Parameters (Sensitivity Variables) ---
candidate_buses = [14,15]; % List of positions to explore
Pgrid_con_values = 900;            % Base size
alpha_values    = [1.0]; % Grid connection ratio Pgrid=alpha*Pel
beta_values     = [1.0]; % Electrolysis ratio Pel=beta*Pwind
cost_c1_wel     = 0;                    % Cost of selling power $/MW
R_wel_val       = 0.7500;                % Droop for WEL
C_RD_wel        = 5;                     % Redispatch cost of WEL $/MW

% Initialize storage for results
results_log = struct();

% --- Pre-Loop Setup (Grid & Static Generators) ---
mpc_base = My_case24_ieee_rts(); 
ng = size(mpc_base.gen,1); nb = size(mpc_base.bus,1); nl = size(mpc_base.branch, 1);

% Generator indexes and droop coefficients
G_f_hydro = [25:30]'; R_hydro = 6*ones(size(G_f_hydro,1),1);
G_f_coal1 = [33]; R_coal1 = 0.8571*ones(size(G_f_coal1,1),1);
G_f_coal2 = [21;22;31;32]; R_coal2 = 1.9355*ones(size(G_f_coal2,1),1);

% Add static onshore wind
P_WPP_MW = 200*ones(6,1); wind_buses = [3;5;7;16;21;23]; cost_c1_wind = 0;
mpc_base = add_wind_gens_from_pu(mpc_base, wind_buses, P_WPP_MW, cost_c1_wind);

%Load duration curve
LOPT= Create_RTS24_LDC();
ldc.P_MW = LOPT(:,1);
ldc.u    = LOPT(:,2);

% Sampler and MC options settings
sampler.method    = 'inv_cdf';
sampler.sigma_rel = 0.04;
sampler.clip_mult = [0.85 1.15];
sampler.clip_m    = [];        
sampler.keep_sum  = true;
sampler.scale_QD  = true;
mpopt = mpoption('model','DC','opf.dc.solver','MIPS','verbose',0);
opts = struct();

% ---- Monte Carlo stopping & accuracy ----
opts.Nmin     = 100;        % minimum valid samples before stopping is allowed
opts.Nmax     = 20000;      % hard cap
opts.rel_tol  = 0.01;       % 1% relative half-width
opts.abs_tol  = 1e-3;       % absolute half-width (CI units)
opts.conf     = 0.95;       % confidence level
opts.print_every = 50;      % progress print frequency

% ---- Congestion index mode ----
opts.mode = 'max';          % 'max' or 'mean' (passed to csopf_droop_CI)

% ---- Wind Monte Carlo settings ----
opts.wind = struct();
opts.wind.enable = true;    % turn wind sampling ON
% KDE / probabilistic models for wind availability
models = build_wind_models_pooled("WindScen.mat");
% Must be a struct array with .sample(n) method
opts.wind.models = models; % <-- output of build_wind_models_pooled()
opts.wind.clip_pu = [0 1];  % enforce physical limits

% enable WEL sampling from FCOPT tables
opts.wel = struct();
opts.wel.enable    = true;
opts.wel.fcopt_file = {'DD_FCOPT_cell.mat'};
opts.wel.clip_pu   = [0 1];
opts.wel.turbine   = {'IEA'};
nwel=numel(opts.wel.fcopt_file);

for k = 1:nwel
    % Get the specific filename for this unit
    this_file = opts.wel.fcopt_file{k};
    % Load the file
    S = load(this_file);
    % Nest the specific table family into the k-th cell
    if isfield(S, 'FCOPT_cell')
        opts.wel.FCOPT_cell{k} = S.FCOPT_cell;
    elseif isfield(S, 'DD_FCOPT_cell')
        opts.wel.FCOPT_cell{k} = S.DD_FCOPT_cell;
    elseif isfield(S, 'DFIG_FCOPT_cell')
        opts.wel.FCOPT_cell{k} = S.DFIG_FCOPT_cell;
    else
        error('Cannot find a valid FCOPT variable in %s', this_file);
    end
end
% ---- (optional) diagnostics ----
opts.debug = false;         % set true if you want extra prints later

% --- Loop over Candidate WEL Positions ---
for i = 1:length(candidate_buses)
    current_bus = candidate_buses(i); % change to i if needed
    Pgrid_con = Pgrid_con_values(1);  % change to i if needed
    alpha_val = alpha_values(1); beta_val = beta_values(1); % change to i if needed

    fprintf('\n--- Testing WEL at Bus %d (%d of %d) ---\n', current_bus, i, length(candidate_buses));
    
    % Copy base case to avoid accumulating generators from previous iterations
    mpc = mpc_base;
    
    % Add the Offshore Wind-Electrolyzer module at the current bus
    % Note: alpha and beta are passed as scalars/vectors based on your function needs
    mpc = add_wel_net_gen(mpc, current_bus, Pgrid_con, cost_c1_wel, alpha_val, beta_val); 
    
    % Update Generator Indices and Droop for this specific WEL location
    G_f_wel = mpc.wind_elec.idx'; 
    R_wel   = R_wel_val * ones(size(G_f_wel,1),1);
    
    mpc.csopf.Gf = [G_f_hydro; G_f_coal1; G_f_coal2; G_f_wel];
    mpc.csopf.R  = [R_hydro; R_coal1; R_coal2; R_wel];
    
    % Define remaining csopf fields (flex_buses, redispatch, costs, etc.)
    mpc.csopf.flex_buses = [];
    mpc.csopf.PLF_max     = [];   
    mpc.csopf.redispatch_gen_idx = [1,2,5,6,21,22,31,32,G_f_wel'];
    mpc.csopf.C_RD = [15,15,15,15,16,16,14,14,C_RD_wel];
    mpc.csopf.C_LS = 1000 * ones(nb,1);
    mpc.csopf.C_LF = 5 * ones(numel(mpc.csopf.flex_buses),1);

    mpc.csopf.Lf     = [];     % bus indices of droop loads
    mpc.csopf.RL     = [];     % droop R for those loads (same length as Lf)
    mpc.csopf.PD_min = [];     % MW min allowed for those loads
    mpc.csopf.PD_max = [];     % MW max allowed  (or equal to PD if no upper)
    
    % Contingency Setup & Ranking
    % (We re-rank per location as the WEL location changes the grid impact)
    gen_prob = [0.1	0.1	0.02 0.02 0.1 0.1 0.02 0.02 0.04 0.04 0.04 0.05	0.05 ... 
            0.05 0 0.02 0.02 0.02	0.02 0.02 0.04 0.04 0.12 0.12 0.01 0.01 ...
        	0.01	0.01	0.01	0.01	0.04	0.04	0.08]; % FOR IEEE 96

    for kg = 1:ng
        mpc.csopf.cont(kg).type = 'inj';
        mpc.csopf.cont(kg).dPd = zeros(nb,1);
        mpc.csopf.cont(kg).dPg_idx = kg;
        mpc.csopf.cont(kg).prob = gen_prob(kg);
    end
    
    % % ---- if including line failures ----
    % branch_dist = [ ...
    % 4.828 88.514 35.406 53.108 80.467 49.890 0.000 43.452 ...
    % 37.015 25.750 25.750 69.202 69.202 0.000 0.000 0.000 ...
    % 0.000 53.108 46.671 53.108 107.826 96.561 43.452 19.312 ...
    % 54.718 54.718 57.936 28.968 25.750 16.093 117.482 28.968 ...
    % 28.968 44.257 44.257 24.140 24.140 75.639]; % cctkm distance of OHL
    
    % OHL_lambda = 0.0022*branch_dist; OHL_mu= 8760/8; % Lines rates in /y
    % branch_prob = OHL_lambda ./ (OHL_lambda + OHL_mu); % OHL unavailabilities
    % 
    % TX_lambda = 0.037; TX_mu=8760/1580; % transformer rates in /y
    % isTX = mpc.branch(:,9) ~= 0;   % tap ratio column
    % branch_prob(isTX) = TX_lambda ./ (TX_lambda + TX_mu); % branch unavailabilities
    
    % k0 = numel(mpc.csopf.cont);   % start index (in case you already added gens)
    
    % for kl = 1:nl
    %     % skip lines already out of service
    %     if mpc.branch(kl, 11) == 0
    %         continue;
    %     end
    % 
    %     k0 = k0 + 1;
    % 
    %     mpc.csopf.cont(k0).type = 'line';
    %     mpc.csopf.cont(k0).line = kl;        % line index (external)
    %     mpc.csopf.cont(k0).prob = branch_prob(kl);       % assign probability
    % end

    [top_idx, ~, ~, ~, ~] = rank_contingencies(mpc, 33); % To select only the most harmful (33 includes all)
    mpc.csopf.cont = mpc.csopf.cont(top_idx);

    % --- Monte Carlo Simulation ---
    % Reset opts for each run if necessary
    [MC_res, sim_results, CI_val] = mc_expected_droop_CI(mpc, ldc, sampler, mpopt, opts);

    % % ---- for debugging ----
    % mpopt = mpoption();
    % results = run_csopf_dc(mpc_s, mpopt); % Uses MIPS for DC
    % CI=csopf_droop_CI(results, mpc_s,'max');
    % csopf_simple_report(results, mpc_s); % Report to debug results
    
    % Store results indexed by the bus number
    results_log(i).bus = current_bus;
    results_log(i).meanE = MC_res.meanE;
    results_log(i).halfCI = MC_res.halfCI;
    results_log(i).n_valid = MC_res.n_valid;
    results_log(i).E_RDp = MC_res.E_RDp;
    results_log(i).E_RDm = MC_res.E_RDm;
    results_log(i).E_LFp = MC_res.E_LFp;
    results_log(i).E_LS = MC_res.E_LS;
    
    fprintf('Result for Bus %d: E[CI]=%.5f\n', current_bus, MC_res.meanE);
end

% --- Final Comparison ---
disp('Simulation Summary:');
disp(struct2table(results_log));
