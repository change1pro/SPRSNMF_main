function [X, Y, info] = load_dataset(rootDir, datasetName)
%LOAD_DATASET Load a .mat dataset and return X as features-by-samples.

dataFile = fullfile(rootDir, [datasetName '.mat']);
if ~exist(dataFile, 'file')
    error('Dataset file not found: %s', dataFile);
end

s = load(dataFile);
if isfield(s, 'X')
    rawX = double(s.X);
else
    error('Dataset %s does not contain variable X.', datasetName);
end

if isfield(s, 'Y')
    Y = double(s.Y(:));
elseif isfield(s, 'y')
    Y = double(s.y(:));
elseif isfield(s, 'gnd')
    Y = double(s.gnd(:));
else
    error('Dataset %s does not contain labels Y/y/gnd.', datasetName);
end

if size(rawX, 1) == numel(Y)
    X = rawX';
else
    X = rawX;
end

X = normalize_data(X);
[~, ~, Y] = unique(Y);

info = struct();
info.name = datasetName;
info.samples = numel(Y);
info.features = size(X, 1);
info.classes = numel(unique(Y));
end
