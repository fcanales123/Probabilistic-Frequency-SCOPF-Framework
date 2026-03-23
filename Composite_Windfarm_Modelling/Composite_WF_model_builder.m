turbine='IEA'; % Turbine model used 'IEA' or 'DTU'
if strcmp(turbine,'IEA')
    WPP_CAP_loop = 0:15:2010;
elseif strcmp(turbine,'DTU')
    WPP_CAP_loop = 0:10:2000;
end
ConvUSD = 1.3; %Conversion ratio pounds to USD
offshore_distance_values = 100; % Renamed to avoid confusion with the scalar in the loop

offshore_site = 'IJ_225'; 

if strcmp(turbine,'IEA')
    num_wind_states_partial=8; % Number of states rounded to in the WT (wind speed) model
elseif strcmp(turbine,'DTU')
    num_wind_states_partial=7; % Number of states rounded to in the WT (wind speed) model
end

num_wind_states_full=14; % Number of states rounded to in the WT (wind speed) model 

num_wind_dir_states = 12; % Number of states of wind direction

num_rounding_states = 'exact'; % Number of states rounded to in the WF model use 'exact' if no rounding required

iteration=0;


% Pre-allocate arrays for efficiency.
% Dimensions are (length(WPP_CAP_loop), length(offshore_distance_values))
EWEA = zeros(length(WPP_CAP_loop), length(offshore_distance_values));
Capacity_factor = zeros(length(WPP_CAP_loop), length(offshore_distance_values));
Distance = zeros(length(WPP_CAP_loop), length(offshore_distance_values));
WPP_OM_COST = zeros(length(WPP_CAP_loop), length(offshore_distance_values));

for x = 1:length(WPP_CAP_loop)
    % Initialize temporary arrays for this parfor iteration
    % These temporary arrays will hold the results for all 'y' values
    % for the current 'x'.
    temp_EWEA = zeros(1, length(offshore_distance_values));
    temp_Capacity_factor = zeros(1, length(offshore_distance_values));
    temp_Distance = zeros(1, length(offshore_distance_values));
    temp_WPP_OM_COST = zeros(1, length(offshore_distance_values));

    for y = 1:length(offshore_distance_values)
        % WPP Parameters
        WPP_CAP = WPP_CAP_loop(x); % Total connected wind farm power in MW
        offshore_distance = offshore_distance_values(y); % Now offshore_distance is a scalar
        
        [rounded_FCOPT, FCOPT_TM_rounded, WT_States,rounded_w_states,rounded_p_states,farm, perBin, diag, WF_States, WF_transition_matrix, FCOPT, FCOPT_TM] = Wind_Farm_Probabilistic_Modelling(turbine, WPP_CAP, offshore_distance, num_wind_states_partial, num_wind_states_full, num_wind_dir_states, num_rounding_states, offshore_site);
        
        temp_EWEA(y) = sum(rounded_FCOPT(:,1) .* rounded_FCOPT(:,2));
        temp_Capacity_factor(y) = temp_EWEA(y) / WPP_CAP;
        temp_Distance(y) = offshore_distance;
        
        % WPP maintenance costs in $/MWh
        if strcmp(turbine,'IEA')
            if offshore_distance < 50 && offshore_distance >= 10
                temp_WPP_OM_COST(y) = (temp_EWEA(y))*(365*24/1000)*ConvUSD*(16.98+0.048*(offshore_distance-10));
            elseif offshore_distance >= 50
                temp_WPP_OM_COST(y) = (temp_EWEA(y))*(365*24/1000)*ConvUSD*(18.9+0.437*(offshore_distance-50));
            end
        elseif strcmp(turbine,'DTU')
            if offshore_distance < 50 && offshore_distance >= 10
                temp_WPP_OM_COST(y) = (temp_EWEA(y)*(365*24/1000)*ConvUSD*(29.53+0.065*(offshore_distance-10)));
            elseif offshore_distance >= 50
                temp_WPP_OM_COST(y) = (temp_EWEA(y)*(365*24/1000)*ConvUSD*(32.13+0.552*(offshore_distance-50)));
            end
        end
    end
    % After the inner loop completes for a given 'x', assign the temporary
    % results to the corresponding slice of the main pre-allocated arrays.
    EWEA(x,:) = temp_EWEA;
    Capacity_factor(x,:) = temp_Capacity_factor;
    Distance(x,:) = temp_Distance;
    WPP_OM_COST(x,:) = temp_WPP_OM_COST;
    % lambda_wT=temp_lambda_wT;
    % mu_wT=temp_mu_wT;
    iteration=iteration+1;
    fprintf('Finished Iteration = %.2f', iteration);
end

%%
ConvUSD = 1.3; %Conversion ratio pounds to USD
turbine='IEA'; % Turbine model used 'IEA' or 'DTU'
WPP_CAP = 2010;
offshore_distance = 100; % Renamed to avoid confusion with the scalar in the loop
num_rounding_states = 5; % Number of states rounded to in the WF model use 'exact' if no rounding required
num_wind_states=5; % Number of states rounded to in the WT (wind speed) model

[rounded_FCOPT, FCOPT_TM_rounded, WT_States] = Wind_Farm_Probabilistic_Modelling(turbine, WPP_CAP, offshore_distance, num_wind_states, num_rounding_states);

