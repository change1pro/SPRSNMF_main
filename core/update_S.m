function S = update_S(F, V, alpha, lambda, gamma, epsv)
%UPDATE_S Row-wise simplex projection for the adaptive similarity matrix.

n = size(F, 1);
if size(V, 1) ~= n
    if size(V, 2) == n
        V = V';
    else
        error('update_S:SizeMismatch', ...
            'F has %d rows, but V is %d-by-%d. V must be samples-by-clusters.', ...
            n, size(V, 1), size(V, 2));
    end
end

distF = pairwise_sqdist(F, F);
distV = pairwise_sqdist(V, V);
if ~isequal(size(distF), size(distV))
    error('update_S:DistanceSizeMismatch', ...
        'distF is %d-by-%d, but distV is %d-by-%d.', ...
        size(distF, 1), size(distF, 2), size(distV, 1), size(distV, 2));
end
Dij = alpha * distF + (lambda / 2) * distV;

if gamma <= epsv
    S = zeros(n, n);
    for ii = 1:n
        row = Dij(ii, :);
        row(ii) = inf;
        [~, jj] = min(row);
        S(ii, jj) = 1;
    end
    return;
end

S = zeros(n, n);
for i = 1:n
    z = -Dij(i, :)' / (2 * gamma);
    z(i) = -inf;
    mask = true(n, 1);
    mask(i) = false;
    si = zeros(n, 1);
    si(mask) = project_simplex(z(mask));
    si(i) = 0;
    S(i, :) = si';
end
S(S < epsv) = 0;
end


