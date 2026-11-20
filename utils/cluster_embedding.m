function pred = cluster_embedding(F, c)
%CLUSTER_EMBEDDING Cluster rows of an embedding matrix by K-means.

data = normalize_rows_local(F);
try
    pred = kmeans(data, c, 'Replicates', 10, 'MaxIter', 1000, 'Display', 'off');
catch
    pred = simple_kmeans_local(data, c, 1000);
end
pred = pred(:);
end

function Xn = normalize_rows_local(X)
rowNorm = sqrt(sum(X.^2, 2));
rowNorm(rowNorm == 0) = 1;
Xn = X ./ rowNorm;
Xn(~isfinite(Xn)) = 0;
end

function pred = simple_kmeans_local(data, c, maxIter)
n = size(data, 1);
idx = randperm(n, c);
centers = data(idx, :);
pred = ones(n, 1);
for t = 1:maxIter
    dist = pairwise_sqdist_local(data, centers);
    [~, newPred] = min(dist, [], 2);
    if isequal(newPred, pred)
        break;
    end
    pred = newPred;
    for j = 1:c
        if any(pred == j)
            centers(j, :) = mean(data(pred == j, :), 1);
        else
            centers(j, :) = data(randi(n), :);
        end
    end
end
end

function D = pairwise_sqdist_local(A, B)
aa = sum(A.^2, 2);
bb = sum(B.^2, 2)';
D = aa + bb - 2 * (A * B');
D = max(D, 0);
end
