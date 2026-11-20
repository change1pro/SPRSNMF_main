function metrics = clustering_metrics(trueLabel, predLabel)
%CLUSTERING_METRICS Compute ACC, NMI, and Purity.

trueLabel = trueLabel(:);
predLabel = predLabel(:);
[~, ~, trueLabel] = unique(trueLabel);
[~, ~, predLabel] = unique(predLabel);

metrics = struct();
metrics.ACC = clustering_acc(trueLabel, predLabel);
metrics.NMI = clustering_nmi(trueLabel, predLabel);
metrics.Purity = clustering_purity(trueLabel, predLabel);
end

function acc = clustering_acc(y, pred)
classes = unique(y);
clusters = unique(pred);
C = zeros(numel(clusters), numel(classes));
for i = 1:numel(clusters)
    for j = 1:numel(classes)
        C(i, j) = sum(pred == clusters(i) & y == classes(j));
    end
end

if exist('matchpairs', 'file') == 2
    try
        pairs = matchpairs(C, 0, 'max');
        total = 0;
        for i = 1:size(pairs, 1)
            total = total + C(pairs(i, 1), pairs(i, 2));
        end
    catch
        total = greedy_match(C);
    end
else
    total = greedy_match(C);
end
acc = total / numel(y);
end

function total = greedy_match(C)
total = 0;
M = C;
while ~isempty(M) && any(M(:) > 0)
    [mx, idx] = max(M(:));
    if mx <= 0
        break;
    end
    [i, j] = ind2sub(size(M), idx);
    total = total + mx;
    M(i, :) = -inf;
    M(:, j) = -inf;
end
end

function nmi = clustering_nmi(y, pred)
n = numel(y);
classes = unique(y);
clusters = unique(pred);
Pij = zeros(numel(classes), numel(clusters));
for i = 1:numel(classes)
    for j = 1:numel(clusters)
        Pij(i, j) = sum(y == classes(i) & pred == clusters(j));
    end
end
Pi = sum(Pij, 2);
Pj = sum(Pij, 1);

MI = 0;
for i = 1:numel(classes)
    for j = 1:numel(clusters)
        if Pij(i, j) > 0
            MI = MI + Pij(i, j) / n * log((n * Pij(i, j)) / (Pi(i) * Pj(j)));
        end
    end
end
Hy = -sum((Pi(Pi > 0) / n) .* log(Pi(Pi > 0) / n));
Hp = -sum((Pj(Pj > 0) / n) .* log(Pj(Pj > 0) / n));
nmi = MI / sqrt(max(Hy * Hp, eps));
end

function purity = clustering_purity(y, pred)
clusters = unique(pred);
total = 0;
for i = 1:numel(clusters)
    idx = pred == clusters(i);
    labels = y(idx);
    total = total + max(histcounts(labels, 1:(max(y) + 1)));
end
purity = total / numel(y);
end
