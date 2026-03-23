function [LOPT] = Create_RTS24_LDC()
    %% Cumulative Load Model Approach
    % This approach typically uses the Load Duration Curve (LDC) for reliability calculations.
    % RTS 24 load parameters cimmulative load model (hourly)
    % Weekly peak load percentages (first 52 weeks)
    Annual_peak=2850;
    weekly_peaks = [ ...
        86.2, 90, 87.8, 83.4, 88, 84.1, 83.2, 80.6, 74, 73.7, 71.5, 72.7, ...
        70.4, 75, 72.1, 80, 75.4, 83.7, 87, 88, 85.6, 81.1, 90, 88.7, 89.6, ...
        86.1, 75.5, 81.6, 80.1, 88, 72.2, 77.6, 80, 72.9, 72.6, 70.5, 78, ...
        69.5, 72.4, 72.4, 74.3, 74.4, 80, 88.1, 88.5, 90.9, 94, 89, 94.2, ...
        97, 100, 95.2];
    
    % Define week categories
    winter_weeks = [1:8, 44:52];
    summer_weeks = 18:30;
    spring_fall_weeks = [9:17, 31:43];
    
    % Hourly profiles for each type (weekday and weekend)
    profiles.winter.weekday = [67 63 60 59 59 60 74 86 95 96 96 95 ...
                               95 95 93 94 99 100 100 96 91 83 73 63];
    profiles.winter.weekend = [78 72 68 66 64 65 66 70 80 88 90 91 ...
                               90 88 87 87 91 100 99 97 94 92 87 81];
    
    profiles.summer.weekday = [64 60 58 56 56 58 64 76 87 95 99 100 ...
                               99 100 100 97 96 96 93 92 92 93 87 72];
    profiles.summer.weekend = [74 70 66 65 64 62 62 66 81 86 91 93 ...
                               93 92 91 91 92 94 95 95 100 93 88 80];
    
    profiles.spring_fall.weekday = [63 62 60 58 59 65 72 85 95 99 100 99 ...
                                     93 92 90 88 90 92 96 98 96 90 80 70];
    profiles.spring_fall.weekend = [75 73 69 66 65 65 68 74 83 89 92 94 ...
                                     91 90 90 86 85 88 92 100 97 95 90 85];
    
    % Preallocate hourly load vector
    hourly_load = zeros(52 * 7 * 24, 1);  % 52 weeks * 7 days * 24 hours
    hour_idx = 1;
    
    for w = 1:52
        % Determine season
        if ismember(w, winter_weeks)
            profile = profiles.winter;
        elseif ismember(w, summer_weeks)
            profile = profiles.summer;
        else
            profile = profiles.spring_fall;
        end
    
        % Get peak for this week
        peak = weekly_peaks(w);
    
        for d = 1:7
            if d <= 5
                daily_profile = profile.weekday;
            else
                daily_profile = profile.weekend;
            end
    
            % Scale profile by weekly peak and add to total load
            hourly_load(hour_idx:hour_idx+23) = peak * daily_profile / 100;
            hour_idx = hour_idx + 24;
        end
    end
    
    % Add day 7 of week 52 (repeat Sunday)
    w = 52;
    peak = weekly_peaks(w);
    
    if ismember(w, winter_weeks)
        profile = profiles.winter;
    elseif ismember(w, summer_weeks)
        profile = profiles.summer;
    else
        profile = profiles.spring_fall;
    end
    
    % Sunday (day 7) is a weekend
    hourly_load(hour_idx:hour_idx+23) = peak * profile.weekend / 100;
    
    % Sort descending to generate load duration curve
    sorted_load = sort(hourly_load, 'descend');
    scaled_load=sorted_load*Annual_peak/100;
    
    % Total hours in the year
    total_hours = 365 * 24;
    total_hours = linspace(0, 1, total_hours);
    
    LOPT=[scaled_load,total_hours.'];
end
