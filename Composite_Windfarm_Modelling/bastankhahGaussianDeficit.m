function [deficit, wakeRadius, uc] = bastankhahGaussianDeficit(D, dw, cw, ct, WS_ref, varargin)
% (same header as before)

% --- params ---
p = inputParser;
addParameter(p, 'k', 0.0324555);
addParameter(p, 'ceps', 0.2);
addParameter(p, 'ctlim', 0.899);
parse(p, varargin{:});
k     = p.Results.k;
ceps  = p.Results.ceps;
ctlim = p.Results.ctlim;

% Ensure arrays
D = D(:);
ct = ct(:);

% 1) Epsilon
ct_clipped = min(ct, ctlim);
sqrt1ct = sqrt(1 - ct_clipped);
beta = 0.5 * (1 + sqrt1ct) ./ sqrt1ct;
epsilon = ceps * sqrt(beta);

% 2) Sigma
sigma = k .* dw + epsilon .* D;
sigma_sqr = sigma.^2;

% 3) Centerline deficit (safe)
sigma_ratio = sigma ./ D;
term = ct ./ (8 .* sigma_ratio.^2);
inside = max(0, 1 - term);   % prevent negatives
deficit_centre = WS_ref .* (1 - sqrt(inside));

% 4) Gaussian profile
exponent = -0.5 .* (cw.^2) ./ sigma_sqr;
deficit = deficit_centre .* exp(exponent);

% 5) Wake radius
wakeRadius = 2 .* sigma;

% 6) Convection velocity
uc = WS_ref .* (1 - deficit_centre ./ WS_ref);
end
