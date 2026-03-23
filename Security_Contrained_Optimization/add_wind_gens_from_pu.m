function mpc = add_wind_gens_from_pu(mpc, wind_buses, P_rated_MW, cost_c1)
%ADD_WIND_GENS_FROM_PU Append wind generators (availability in pu) to MATPOWER case.
%
% mpc = add_wind_gens_from_pu(mpc, wind_buses, pw_pu, P_rated_MW, cost_c1)
%
% Inputs
%   wind_buses : Nx1 bus indices (external numbering) where wind plants connect
%   pw_pu      : Nx1 availability in pu (0..1 typically) for each plant
%   P_rated_MW : scalar rated power per plant (default 200)
%   cost_c1    : optional linear $/MW cost (default 0). Use small value if desired.
%
% Behavior
%   Pg = P_avail, Pmax = P_avail, Pmin = 0, status = 1
%   Q limits set to 0 (DC OPF ignores Q anyway)

if nargin < 4 || isempty(P_rated_MW), P_rated_MW = 200; end
if nargin < 5 || isempty(cost_c1),    cost_c1 = 0;      end

define_constants;

N = numel(wind_buses);

% ---- ensure bus indices are valid ----
nb = size(mpc.bus, 1);
if any(wind_buses < 1) || any(wind_buses > nb)
    error('Some wind_buses are out of range 1..%d (external bus order).', nb);
end

% ---- available power in MW (clip to [0, P_rated]) ----
P_avail_MW = 0; % Let the wind generators not participate in ranking of contingencies

% ---- build gen rows ----
% gen format:
% [BUS PG QG QMAX QMIN VG MBASE STATUS PMAX PMIN ... ramp/apf]
ng_old = size(mpc.gen, 1);
gen_new = zeros(N, size(mpc.gen, 2));

gen_new(:, GEN_BUS)     = wind_buses;
gen_new(:, PG)          = P_avail_MW;     % initial Pg
gen_new(:, QG)          = 0;
gen_new(:, QMAX)        = 0;
gen_new(:, QMIN)        = 0;
gen_new(:, VG)          = 1.0;
gen_new(:, MBASE)       = mpc.baseMVA;    % can also use 100
gen_new(:, GEN_STATUS)  = 1;
gen_new(:, PMAX)        = P_avail_MW;     % cannot exceed available wind
gen_new(:, PMIN)        = 0;              % allow curtailment

% Optional: if gen matrix has extra columns, leave them zeros (OK)

% ---- append generators ----
mpc.gen = [mpc.gen; gen_new];

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
if ~isfield(mpc, 'wind')
    mpc.wind = struct();
end
mpc.wind.idx = (ng_old+1):(ng_old+N);
mpc.wind.buses = wind_buses;
mpc.wind.P_rated_MW = P_rated_MW;

end
