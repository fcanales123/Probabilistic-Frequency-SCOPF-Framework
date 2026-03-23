function models = build_wind_models_pooled(matfile, opts)
%BUILD_WIND_MODELS_POOLED  One PDF per plant using all hours pooled.
%
% models(p) fields:
%   .support    [lo hi]
%   .x          grid points
%   .pdf        pdf(x)
%   .cdf        cdf(x) (normalized)
%   .inv_cdf    function handle u->x via inverse CDF
%   .sample     function handle n->samples (inverse transform)
%   .n_samples  number of data points used
%
% Expected input WindScen.mat contains either:
%   - WindScen : 1xNplants cell, each cell is [T x Ns] (e.g., 43 x 100)
%     OR
%   - a numeric array [T x Ns x Nplants] or [Nplants x T x Ns] (we try to detect)
%
% opts (optional):
%   opts.n_grid   = 400 (grid size)
%   opts.support  = []  (auto) or [lo hi]
%   opts.bandwidth = [] (let ksdensity choose) or numeric
%   opts.pad_frac = 1e-6 (support padding fraction)
%
% Example:
%   models = build_wind_models_pooled("WindScen.mat");
%   x = models(3).sample(10000);

if nargin < 2, opts = struct(); end
if ~isfield(opts, 'n_grid'),     opts.n_grid = 400; end
if ~isfield(opts, 'support'),    opts.support = []; end
if ~isfield(opts, 'bandwidth'),  opts.bandwidth = []; end
if ~isfield(opts, 'pad_frac'),   opts.pad_frac = 1e-6; end

S = load(matfile);

% --- load dataset into a standard cell array: WindScen{p} = [T x Ns] ---
if isfield(S, 'WindScen') && iscell(S.WindScen)
    WindScen = S.WindScen;
elseif isfield(S, 'WindScen')
    X = S.WindScen;
    WindScen = numeric_to_cell_plants(X);
else
    % fallback: pick the first variable in the file
    fn = fieldnames(S);
    X = S.(fn{1});
    if iscell(X)
        WindScen = X;
    else
        WindScen = numeric_to_cell_plants(X);
    end
end

nP = numel(WindScen);
models = repmat(struct(), nP, 1);

for p = 1:nP
    Y = WindScen{p};                % [T x Ns]
    y = Y(:);
    y = y(isfinite(y));             % remove NaN/Inf

    if isempty(y)
        error('Plant %d: no finite data points found.', p);
    end

    % --- support ---
    if isempty(opts.support)
        lo = min(y);
        hi = max(y);
        if lo == hi
            lo = lo - 1e-6; hi = hi + 1e-6;
        end
        pad = opts.pad_frac * max(1, abs(hi - lo));
        sup = [lo - pad, hi + pad];
    else
        sup = opts.support(:).';
        if numel(sup) ~= 2 || sup(1) >= sup(2)
            error('opts.support must be [low high] with low < high, or [].');
        end
    end

    % Clip data to support (helps ksdensity avoid support errors)
    y = min(max(y, sup(1)), sup(2));

    % --- grid ---
    x = linspace(sup(1), sup(2), opts.n_grid).';

    % --- KDE ---
    if isempty(opts.bandwidth)
        [pdfx, xout] = ksdensity(y, x, 'Support', sup, 'Function', 'pdf');
    else
        [pdfx, xout] = ksdensity(y, x, 'Support', sup, 'Function', 'pdf', ...
            'Bandwidth', opts.bandwidth);
    end

    % --- CDF (numerical integration + normalize) ---
    cdfx = cumtrapz(xout, pdfx);
    if cdfx(end) <= 0 || ~isfinite(cdfx(end))
        error('Plant %d: invalid CDF normalization.', p);
    end
    cdfx = cdfx / cdfx(end);
    cdfx(1) = 0; cdfx(end) = 1;

    % --- make CDF strictly usable for inverse transform ---
    % enforce monotone nondecreasing (guards tiny numerical wiggles)
    cdfx = max(cdfx, [0; cdfx(1:end-1)]);
    
    % keep only unique CDF points for interp1
    [cdfu, ia] = unique(cdfx, 'stable');
    xu = xout(ia);
    
    % ensure endpoints exist
    if cdfu(1) > 0
        cdfu = [0; cdfu];
        xu   = [xu(1); xu];
    end
    if cdfu(end) < 1
        cdfu = [cdfu; 1];
        xu   = [xu; xu(end)];
    end
    
    inv_cdf = @(u) interp1(cdfu, xu, u, 'linear', 'extrap');
    sampler = @(n) min(max(inv_cdf(rand(n,1)), sup(1)), sup(2));

    % store
    models(p).plant     = p;
    models(p).support   = sup;
    models(p).x         = xout;
    models(p).pdf       = pdfx;
    models(p).cdf       = cdfx;
    models(p).inv_cdf   = inv_cdf;
    models(p).sample    = sampler;
    models(p).n_samples = numel(y);
end
end

% ---------- helper ----------
function C = numeric_to_cell_plants(X)
% Try to convert numeric array to {plant} cells of [T x Ns]
sz = size(X);

% Common: [T x Ns x Nplants]
if numel(sz) == 3
    if sz(3) >= 2
        T = sz(1); Ns = sz(2); nP = sz(3);
        C = cell(1, nP);
        for p = 1:nP
            C{p} = reshape(X(:,:,p), T, Ns);
        end
        return;
    end
end

% Common: [Nplants x T x Ns]
if numel(sz) == 3 && sz(1) >= 2
    nP = sz(1); T = sz(2); Ns = sz(3);
    C = cell(1, nP);
    for p = 1:nP
        C{p} = reshape(squeeze(X(p,:,:)), T, Ns);
    end
    return;
end

error('Cannot interpret numeric WindScen array dimensions: %s', mat2str(sz));
end
