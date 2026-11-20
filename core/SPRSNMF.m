function [F, obj] = SPRSNMF(X, c, alpha, lambda, beta, k)
%SPRSNMF Self-paced adaptive graph sparse NMF.

%% 初始化
maxIter = 80;
innerB = 5;
innerF = 5;
tol = 1e-5;
epsv = 1e-10;
initialPercent = 50;
percentStep = 10;
zetaRatio = 0.5;
initNMFIter = 20;

[~, n] = size(X);
label = kmeans(X', c, 'Replicates', 10, 'MaxIter', 1000, 'Display', 'off');
F = epsv * ones(n, c);
for i = 1:n
    F(i, label(i)) = 1;
end
B = max(X * F / max(F' * F, epsv), epsv);

for t = 1:max(5, round(initNMFIter / 2))
    B = B .* ((X * F) ./ max(B * (F' * F), epsv));
    F = F .* ((X' * B) ./ max(F * (B' * B), epsv));
    B = max(B, epsv);
    F = max(F, epsv);
end
F = F ./ max(sqrt(sum(F.^2, 1)), epsv);

[~, labelF] = max(F, [], 2);
F0 = epsv * ones(n, c);
for i = 1:n
    F0(i, labelF(i)) = max(F(i, labelF(i)), epsv);
end
F = F0 ./ max(sqrt(sum(F0.^2, 1)), epsv);

dist = pairwise_sqdist(F, F);
S = zeros(n, n);
for i = 1:n
    [~, idx] = sort(dist(i, :), 'ascend');
    idx(idx == i) = [];
    idx = idx(1:min(k, n - 1));
    S(i, idx) = 1 / numel(idx);
end

obj = [];

%% Main iteration
for iter = 1:maxIter
    E = X - B * F';
    loss = sqrt(sum(E.^2, 1))';

    pct = min(initialPercent + (iter - 1) * percentStep, 100);
    lossSort = sort(loss);
    eta = max(lossSort(max(1, ceil(pct / 100 * n))), epsv);
    zeta = zetaRatio * eta;

    w = zeros(n, 1);
    low = zeta * eta / (zeta + eta);
    w(loss <= low) = 1;
    mid = loss > low & loss < eta;
    w(mid) = zeta ./ max(loss(mid), eps) - zeta / eta;
    w = min(max(w, 0), 1);
    D = w ./ (2 * sqrt(loss.^2 + epsv));

    W = (S + S') / 2;
    Ds = diag(sum(W, 2));
    Ls = Ds - W;
    [U, Dv] = eig(full((Ls + Ls') / 2));
    [~, idx] = sort(diag(Dv), 'ascend');
    V = U(:, idx(1:c));

    distF = pairwise_sqdist(F, F);
    distV = pairwise_sqdist(V, V);
    Dij = alpha * distF + (lambda / 2) * distV;
    gammaVals = zeros(n, 1);
    kk = min(max(round(k), 1), max(n - 2, 1));
    for i = 1:n
        row = Dij(i, :);
        row(i) = inf;
        row = sort(row, 'ascend');
        gammaVals(i) = kk / 2 * row(kk + 1) - 0.5 * sum(row(1:kk));
    end
    gamma = max(mean(gammaVals), epsv);

    S = update_S(F, V, alpha, lambda, gamma, epsv);
    W = (S + S') / 2;
    Ds = diag(sum(W, 2));

    B = update_B(X, F, D, B, innerB, epsv);
    F = update_f_SPRSNMF(X, B, W, Ds, D, F, alpha, beta, innerF, epsv);

    [~, labelF] = max(F, [], 2);
    F0 = epsv * ones(n, c);
    for i = 1:n
        F0(i, labelF(i)) = max(F(i, labelF(i)), epsv);
    end
    F = F0 ./ max(sqrt(sum(F0.^2, 1)), epsv);

    Q = diag(1 ./ (4 * (sum(F.^2, 2) + epsv).^(3/4)));
    Ls = Ds - W;
    recon = sum(w .* sqrt(sum((X - B * F').^2, 1))');
    graphF = alpha * sum(sum(pairwise_sqdist(F, F) .* S));
    smoothS = gamma * sum(sum(S.^2));
    spec = lambda * trace(V' * Ls * V);
    sparse = beta * trace(F' * Q * F);
    selfPaced = -sum(zeta * log(w + zeta / eta));

    obj(iter) = recon + graphF + smoothS + spec + sparse + selfPaced;

    if iter > 2
        relchg = abs(obj(iter - 1) - obj(iter)) / max(1, abs(obj(iter - 1)));
        if relchg < tol
            break;
        end
    end
end

end
