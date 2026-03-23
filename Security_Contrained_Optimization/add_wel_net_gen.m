function mpc = add_wel_net_gen(mpc, bus_wel, Pgrid_con_MW, cost_c1, alpha, beta)
    define_constants;
    N = numel(bus_wel);
    ng_old = size(mpc.gen, 1);

    % --- build a new gen row (DC OPF only really uses Pg, Pmax, Pmin, status, bus) ---
    newgen = zeros(N, size(mpc.gen, 2));

    newgen(:,GEN_BUS) = bus_wel;
    newgen(:,PG)      = 0;          % initial
    newgen(:,QG)      = 0;
    newgen(:,QMAX)    = 0;
    newgen(:,QMIN)    = 0;
    newgen(:,VG)      = 1;
    newgen(:,MBASE)   = mpc.baseMVA;
    newgen(:,GEN_STATUS) = 1;
    newgen(:,PMAX)    = 0;   % upper export limit (will be set per-sample)
    newgen(:,PMIN)    = 0;

    mpc.gen = [mpc.gen; newgen];

    % ---- append gencost rows ----
    % gencost format (polynomial):
    % [MODEL STARTUP SHUTDOWN NCOST c(n-1) ... c0]
    % For linear cost: NCOST=2 => [c1 c0]
    gc_new = zeros(N, size(mpc.gencost, 2));
    gc_new(:, MODEL)    = POLYNOMIAL;
    gc_new(:, STARTUP)  = 0;
    gc_new(:, SHUTDOWN) = 0;
    gc_new(:, NCOST)    = 2;
    
    % c1 and c0 in the last two columns (MATPOWER convention)
    % If your gencost has more columns, we place them at the end.
    gc_new(:, end-1) = cost_c1;   % $/MW
    gc_new(:, end)   = 0;         % $
    
    mpc.gencost = [mpc.gencost; gc_new];

    % Optional: tag the units (not a MATPOWER field; just info)
    if ~isfield(mpc, 'wind_elec')
        mpc.wind_elec = struct();
    end
    mpc.wind_elec.idx           = (ng_old+1):(ng_old+N);
    mpc.wind_elec.buses         = bus_wel(:);
    mpc.wind_elec.grid_con_MW   = Pgrid_con_MW(:);
    mpc.wind_elec.alpha         = alpha(:);          % N x 1 (or scalar)
    mpc.wind_elec.beta          = beta(:);           % N x 1 (or scalar)
    mpc.wind_elec.Pel_cap_MW    = Pgrid_con_MW./alpha;           % N x 1 (or scalar)
    mpc.wind_elec.Pw_cap_MW     = Pgrid_con_MW./(alpha.*beta);   % N x 1 (or scalar)
    mpc.wind_elec.droop_up      = zeros(N,1);           % N x 1 (or scalar)
    mpc.wind_elec.droop_dwn     = zeros(N,1);           % N x 1 (or scalar)
end
