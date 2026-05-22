function [Pw_avail_pu, Pw_avail_MW, idx] = sample_wel_from_FCOPT(FCOPT_cell, Pw_cap_MW, turbine)
%SAMPLE_WEL_FROM_FCOPT Sample available wind power from the FCOPT table family.
%
%   [Pw_avail_pu, Pw_avail_MW, idx] = sample_wel_from_FCOPT(FCOPT_cell, Pw_cap_MW)
%
% Assumptions:
%   - FCOPT_cell is a 1 x N cell array where N corresponds to capacities 0:15:2010 MW
%   - Pw_cap_MW must be exactly on the grid (0:15:2010). This function enforces it.
%   - FCOPT_cell{1,idx} contains either:
%       (A) numeric matrix [x, p] where p are probabilities (pmf) or weights
%       (B) numeric matrix [x, cdf] where cdf is nondecreasing and ends near 1
%       (C) struct with fields like .x and .p (or .cdf)
%
% Output:
%   Pw_avail_pu in [0,1] typically
%   Pw_avail_MW = Pw_avail_pu * Pw_cap_MW

    % ----- map capacity to table index -----
    if nargin < 3 || isempty(turbine), turbine = 'IEA'; end
    turbine = upper(string(turbine));
    
    switch turbine
        case "IEA"
            step = 15;
            maxCap = 2010;
        case "DTU"
            step = 10;
            maxCap = 2000;
        otherwise
            error('Unknown turbine "%s". Use "IEA" or "DTU".', turbine);
    end
    
    if Pw_cap_MW < 0 || Pw_cap_MW > maxCap
        error('Pw_cap_MW must be within [0, %d] for %s. Got %.4g.', maxCap, turbine, Pw_cap_MW);
    end

    idx = round(Pw_cap_MW/step) + 1;

    if idx < 1 || idx > numel(FCOPT_cell)
        error('Index idx=%d out of bounds for FCOPT_cell (numel=%d).', idx, numel(FCOPT_cell));
    end

    T = FCOPT_cell{1, idx};
    if isempty(T)
        error('FCOPT_cell{1,%d} is empty for Pw_cap_MW=%.4g.', idx, Pw_cap_MW);
    end

    % ----- extract x and probability description -----
    if isstruct(T)
        if isfield(T,'x')
            x = T.x(:);
        elseif isfield(T,'X')
            x = T.X(:);
        else
            error('Struct table at idx=%d has no x/X field.', idx);
        end

        if isfield(T,'p')
            p = T.p(:);
            mode = 'pmf';
        elseif isfield(T,'pdf')
            p = T.pdf(:);
            mode = 'pmf';
        elseif isfield(T,'cdf')
            cdf = T.cdf(:);
            mode = 'cdf';
        else
            error('Struct table at idx=%d has no p/pdf/cdf field.', idx);
        end

    elseif isnumeric(T)
        if size(T,2) < 2
            error('Numeric table at idx=%d must have at least 2 columns [x, p/cdf].', idx);
        end
        x = T(:,1);
        y = T(:,2);

        % decide whether 2nd col looks like cdf
        y = y(:);
        if all(isfinite(y)) && all(diff(y) >= -1e-10) && y(1) >= -1e-6 && y(end) <= 1+1e-6
            % looks like a CDF
            cdf = min(max(y,0),1);
            mode = 'cdf';
        else
            % treat as pmf/weights
            p = y;
            mode = 'pmf';
        end
    else
        error('Unsupported table type at idx=%d: %s', idx, class(T));
    end

    % ----- build CDF and sample -----
    switch mode
        case 'pmf'
            p = p(:);
            x = x(:);

            % clean invalid entries
            mask = isfinite(x) & isfinite(p) & (p >= 0);
            x = x(mask);
            p = p(mask);

            if isempty(x)
                error('After cleaning, table at idx=%d is empty.', idx);
            end

            s = sum(p);
            if s <= 0
                error('Probabilities/weights sum to <=0 at idx=%d.', idx);
            end
            p = p / s;

            cdf = cumsum(p);
            cdf(end) = 1;  % enforce exact

            u = rand();
            j = find(cdf >= u, 1, 'first');
            if isempty(j), j = numel(x); end
            Pw_avail_MW = x(j);

        case 'cdf'
            x = x(:);
            cdf = cdf(:);

            % clean and enforce monotonic
            mask = isfinite(x) & isfinite(cdf);
            x = x(mask);
            cdf = cdf(mask);

            if isempty(x)
                error('After cleaning, CDF table at idx=%d is empty.', idx);
            end

            % make cdf monotone nondecreasing and clipped
            cdf = min(max(cdf,0),1);
            cdf = cummax(cdf);
            if cdf(end) < 1
                % normalize if ends below 1
                if cdf(end) <= 0
                    error('CDF ends at 0 at idx=%d; invalid.', idx);
                end
                cdf = cdf / cdf(end);
            end
            cdf(end) = 1;

            % IMPORTANT: interp1 requires unique sample points
            % enforce uniqueness on cdf by keeping first occurrence
            [cdf_u, ia] = unique(cdf, 'stable');
            x_u = x(ia);

            u = rand();
            Pw_avail_MW = interp1(cdf_u, x_u, u, 'linear', 'extrap');

        otherwise
            error('Internal error: unknown mode "%s".', mode);
    end

    % ----- guard + MW conversion -----
    if ~isfinite(Pw_avail_MW)
        error('Sampled Pw_avail_pu is not finite at idx=%d.', idx);
    end

    % If the table is in MW already, you can modify here.
    % For now we assume x is PU in [0,1].
    Pw_avail_MW = max(Pw_avail_MW, 0);
    Pw_avail_pu = Pw_avail_MW / Pw_cap_MW;
end
