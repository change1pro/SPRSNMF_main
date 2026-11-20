function X = normalize_data(X)
%NORMALIZE_DATA Make data nonnegative and normalize each sample.

X = double(X);
minVal = min(X(:));
if minVal < 0
    X = X - minVal;
end
maxVal = max(X(:));
if maxVal > 1
    X = X / maxVal;
end
colNorm = sqrt(sum(X.^2, 1));
colNorm(colNorm == 0) = 1;
X = X ./ colNorm;
X = max(X, 0);
end
