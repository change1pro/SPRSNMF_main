clear
clc

rootDir = fileparts(mfilename('fullpath'));

addpath(fullfile(rootDir, 'core'), '-begin');
addpath(fullfile(rootDir, 'utils'), '-begin');

datasets = {'JAFFE', 'GLIOMA', 'lung'};

outDir = fullfile(rootDir, 'visualization_results', 'clustering_visualization');
if exist(outDir, 'dir') == 0
    mkdir(outDir);
end

allRows = {};

figure('Color', 'w', 'Position', [80, 40, 960, 980]);
layout = tiledlayout(numel(datasets), 2, ...
    'TileSpacing', 'compact', 'Padding', 'loose');

for d = 1:numel(datasets)
    dataset = datasets{d};
    fprintf('\nDataset: %s\n', dataset);

    [X, Y, info] = load_dataset(rootDir, dataset);
    c = info.classes;
    Y = Y(:);

    [alpha, lambda, beta, k] = read_splagnmf_best_params(rootDir, dataset);
    fprintf('Use parameters: alpha=%g, lambda=%g, beta=%g, k=%g\n', ...
        alpha, lambda, beta, k);

    % Fixed seeds make SPLAGNMF initialization and t-SNE reproducible.
    rng(20260813 + d * 100, 'twister');
    [F, ~] = splagnmf(X, c, alpha, lambda, beta, k);
    pred = cluster_embedding(F, c);

    raw2d = reduce_to_2d(X');
    f2d = reduce_to_2d(F);

    nexttile(layout, 2 * d - 1);
    draw_scatter(raw2d, Y, c);
    title(sprintf('(%c)  the original data %s', char('a' + 2 * d - 2), dataset), ...
        'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');

    nexttile(layout, 2 * d);
    draw_scatter(f2d, pred, c);
    title(sprintf('(%c) the representation of %s obtained by SPRSNMF ', char('a' + 2 * d - 1), dataset), ...
        'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');

    for i = 1:size(raw2d, 1)
        allRows(end + 1, :) = {dataset, i, Y(i), pred(i), ...
            raw2d(i, 1), raw2d(i, 2), f2d(i, 1), f2d(i, 2), ...
            alpha, lambda, beta, k};
    end
end

drawnow;
plotAxes = findall(gcf, 'Type', 'axes');
for a = 1:numel(plotAxes)
    ax = plotAxes(a);
    ax.XTick = ax.XTick;
    ax.YTick = ax.YTick;
    ax.XTickLabel = strrep(ax.XTickLabel, '-', char(8722));
    ax.YTickLabel = strrep(ax.YTickLabel, '-', char(8722));
end

exportgraphics(gcf, fullfile(outDir, 'SPLAGNMF_clustering_visualization.png'), ...
    'Resolution', 300);
savefig(gcf, fullfile(outDir, 'SPLAGNMF_clustering_visualization.fig'));
exportgraphics(gcf, fullfile(outDir, 'SPLAGNMF_clustering_visualization.pdf'), ...
    'ContentType', 'vector');

coordTable = cell2table(allRows, 'VariableNames', ...
    {'Dataset', 'SampleIndex', 'TrueLabel', 'PredLabel', ...
    'RawDim1', 'RawDim2', 'FDim1', 'FDim2', ...
    'alpha', 'lambda', 'beta', 'k'});
writetable(coordTable, fullfile(outDir, 'SPLAGNMF_clustering_visualization_coordinates.csv'));

fprintf('\nClustering visualization saved to:\n%s\n', outDir);

function [alpha, lambda, beta, k] = read_splagnmf_best_params(rootDir, dataset)

file = fullfile(rootDir, 'SPLAGNMF_top_params', [dataset '_top50_params.csv']);

T = readtable(file);
[~, idx] = max(T.Score_mean);

alpha = T.alpha(idx);
lambda = T.lambda(idx);
beta = T.beta(idx);
k = T.k(idx);
end

function Y2 = reduce_to_2d(data)
data = double(data);
data(~isfinite(data)) = 0;

perplexity = min(30, max(5, floor((size(data, 1) - 1) / 3)));
Y2 = tsne(data, 'NumDimensions', 2, 'Perplexity', perplexity, ...
    'Standardize', true);
end

function draw_scatter(Y2, label, c)
[~, ~, label] = unique(label(:));
colors = lines(c);
hold on
for i = 1:max(label)
    idx = label == i;
    scatter(Y2(idx, 1), Y2(idx, 2), 20, colors(i, :), ...
        'filled', 'MarkerFaceAlpha', 0.78, 'MarkerEdgeColor', 'none');
end
hold off
axis tight
box on
grid on
set(gca, 'FontName', 'Times New Roman', ...
    'FontSize', 10, 'FontWeight', 'bold', 'LineWidth', 0.9);
xlabel('Dim 1', 'FontName', 'Times New Roman', ...
    'FontSize', 11, 'FontWeight', 'bold');
ylabel('Dim 2', 'FontName', 'Times New Roman', ...
    'FontSize', 11, 'FontWeight', 'bold');
end
