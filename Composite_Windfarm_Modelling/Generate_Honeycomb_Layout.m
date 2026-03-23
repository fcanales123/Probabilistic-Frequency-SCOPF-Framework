function [X_coords, Y_coords] = Generate_Honeycomb_Layout(nT, turbine, wT_spacing)
%Generate_Honeycomb_Layout Honeycomb (hex) layout with near-square footprint.
%   Balances columns and rows so the farm looks as square as possible while
%   keeping proper hex spacing (staggered columns).

    % --- 1) Rotor diameter by turbine type ---
    switch upper(turbine)
        case 'IEA'
            D = 240;      % IEA 15 MW
        case 'DTU'
            D = 178.3;    % DTU 10 MW
        otherwise
            error('Invalid turbine type. Choose ''IEA'' or ''DTU''.');
    end

    % --- 2) Hex spacing (center-to-center) ---
    S = wT_spacing * D;     % e.g., 7D
    col_dist = S * sqrt(3)/2; % horizontal distance between columns
    row_dist = S;             % vertical distance within a column

    % --- 3) Choose number of columns close to square, distribute rows ---
    n_cols = max(1, round(sqrt(nT)));       % aim for square
    base_rows = floor(nT / n_cols);
    r_extra   = mod(nT, n_cols);            % first r_extra columns get +1 row

    % If base_rows is zero (very small nT), reduce columns
    while base_rows == 0
        n_cols = n_cols - 1;
        base_rows = floor(nT / n_cols);
        r_extra   = mod(nT, n_cols);
    end

    % --- 4) Generate coordinates column by column (staggering every other) ---
    X = [];
    Y = [];
    % X.reserve = []; %#ok<STRNU>  % no-op in MATLAB, here for clarity

    % Center columns around zero
    x0 = -(n_cols-1)/2 * col_dist;

    for c = 1:n_cols
        % Rows in this column: either base_rows or base_rows+1
        rows_c = base_rows + (c <= r_extra);

        % y positions for this column, centered
        y_col = ((0:rows_c-1) - (rows_c-1)/2) * row_dist;

        % Stagger every second column for the honeycomb
        % (You can flip even/odd if you want the opposite staggering.)
        if mod(c,2) == 0
            y_col = y_col + row_dist/2;
        end

        x_c = x0 + (c-1) * col_dist;
        X = [X; repmat(x_c, rows_c, 1)]; %#ok<AGROW>
        Y = [Y; y_col(:)];               %#ok<AGROW>
    end

    % We might have placed exactly nT. If a rounding edge case over-placed, trim.
    if numel(X) > nT
        X = X(1:nT);
        Y = Y(1:nT);
    end

    % --- 5) Center the final layout (should already be centered, just ensure) ---
    X_coords = X - mean(X);
    Y_coords = Y - mean(Y);

    % % --- 6) Quick visualization ---
    % figure;
    % plot(X_coords, Y_coords, 'o', 'MarkerSize', 8, ...
    %     'MarkerFaceColor', [0.1 0.5 0.8], 'MarkerEdgeColor', 'k');
    % hold on;
    % try
    %     viscircles([X_coords(1), Y_coords(1)], D/2, 'Color', 'r', ...
    %         'LineWidth', 1, 'LineStyle', '--');
    % catch
    %     % viscircles requires Image Processing Toolbox; ignore if unavailable
    % end
    % axis equal; grid on;
    % xlabel('X (m)'); ylabel('Y (m)');
    % title(sprintf('Honeycomb Layout: %d %s Turbines (Spacing %.1f m)', nT, upper(turbine), S));
end
