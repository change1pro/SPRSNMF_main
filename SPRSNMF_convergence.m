clear
clc

rootDir = 'D:\matlab-lhb\SPRSNMF';
addpath(fullfile(rootDir, 'utils'), '-begin');
addpath(fullfile(rootDir, 'core'), '-begin');

datasets = {'att','UMIST','JAFFE','oh15','GLIOMA','BA', ...
    'lymphoma','lung','colon'};

coarseDir = fullfile(rootDir, 'SPRSNMF_results_coarse');
outDir = fullfile(rootDir, 'plot_objimg');
if exist(outDir, 'dir') == 0
    mkdir(outDir);
end

figure('Color', 'w', 'Position', [100, 80, 980, 760]);
tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

letters = 'abcdefghijklmnopqrstuvwxyz';
summaryRows = {};

for d = 1:numel(datasets)
    dataset = datasets{d};

    modelFiles = dir(fullfile(coarseDir, dataset, ...
        'alpha-*', 'beta-*', 'lambda-*', 'k-*', 'model.mat'));

    bestScore = -inf;
    bestCurve = [];
    bestFile = '';
    bestInfo = struct('alpha', NaN, 'beta', NaN, 'lambda', NaN, 'k', NaN, ...
        'decreaseRatio', NaN, 'violationRatio', NaN, 'smoothness', NaN);

    for i = 1:numel(modelFiles)
        modelFile = fullfile(modelFiles(i).folder, modelFiles(i).name);
        try
            S = load(modelFile, 'model');
            curve = average_obj_curve(S.model.obj);
        catch
            continue;
        end

        if numel(curve) < 4
            continue;
        end

        [score, info] = convergence_shape_score(curve);
        if score > bestScore
            bestScore = score;
            bestCurve = curve;
            bestFile = modelFile;
            bestInfo = info;
        end
    end

    nexttile;
    if isempty(bestCurve)
        title(sprintf('(%c) %s', letters(d), dataset), ...
            'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');
        text(0.5, 0.5, 'No curve found', 'HorizontalAlignment', 'center');
        axis off
        continue;
    end

    orange = [0.9000, 0.1800, 0.0000];
    plot(1:numel(bestCurve), bestCurve, '-', ...
        'Color', orange, ...
        'LineWidth', 1.8);
    title(sprintf('(%c) %s', letters(d), dataset), ...
        'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Iteration number', 'FontName', 'Times New Roman', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Objective function value', 'FontName', 'Times New Roman', 'FontSize', 11, 'FontWeight', 'bold');
    set(gca, 'FontName', 'Times New Roman', ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'LineWidth', 0.9, ...
        'GridColor', [0.85, 0.85, 0.85], ...
        'GridAlpha', 0.35);
    grid on
    box on
    xlim([1, max(1, numel(bestCurve))]);

    [alphaVal, betaVal, lambdaVal, kVal] = parse_params_from_path(bestFile);
    summaryRows(end + 1, :) = {dataset, bestFile, numel(bestCurve), ...
        bestCurve(1), bestCurve(end), alphaVal, betaVal, lambdaVal, kVal, ...
        bestInfo.decreaseRatio, bestInfo.violationRatio, bestInfo.smoothness, bestScore};

end

pngFile = fullfile(outDir, 'SPRSNMF_convergence_9datasets_bestshape.png');
pdfFile = fullfile(outDir, 'SPRSNMF_convergence_9datasets_bestshape_solid.pdf');
figFile = fullfile(outDir, 'SPRSNMF_convergence_9datasets_bestshape_solid.fig');
csvOut = fullfile(outDir, 'SPRSNMF_convergence_9datasets_bestshape_solid_summary.csv');

exportgraphics(gcf, pngFile, 'Resolution', 300);
exportgraphics(gcf, pdfFile, 'ContentType', 'vector');
savefig(gcf, figFile);

summary = cell2table(summaryRows, 'VariableNames', ...
    {'Dataset', 'ModelFile', 'Iterations', 'StartObj', 'EndObj', ...
    'alpha', 'beta', 'lambda', 'k', ...
    'DecreaseRatio', 'ViolationRatio', 'Smoothness', 'ShapeScore'});
writetable(summary, csvOut);

fprintf('\nSaved best-shaped convergence figure to:\n%s\n%s\n%s\n', pngFile, pdfFile, figFile);

function [score, info] = convergence_shape_score(curve)
curve = curve(:);
curve = curve(isfinite(curve));

span = max(abs(curve(1) - curve(end)), eps);
scale = max(abs(curve(1)), eps);
diffCurve = diff(curve);

decreaseRatio = (curve(1) - curve(end)) / scale;
violationRatio = sum(max(diffCurve, 0)) / span;

if numel(diffCurve) >= 2
    smoothness = sum(abs(diff(diffCurve))) / max(sum(abs(diffCurve)), eps);
else
    smoothness = 0;
end

tailLen = min(5, numel(curve) - 1);
tailChange = mean(abs(diffCurve(end - tailLen + 1:end))) / scale;

    negativePenalty = 0;
    if min(curve) < 0
        negativePenalty = 1000 + abs(min(curve)) / scale;
    end

    weakDecreasePenalty = 0;
    if decreaseRatio <= 0
        weakDecreasePenalty = 1000 + abs(decreaseRatio);
    end

    % A good representative convergence curve should drop clearly, have few
    % upward jumps, and become stable near the end.
    score = 20 * decreaseRatio ...
    - 80 * violationRatio ...
    - 8 * smoothness ...
    - 40 * tailChange ...
    - 100 * negativePenalty ...
    - 100 * weakDecreasePenalty;

info.decreaseRatio = decreaseRatio;
info.violationRatio = violationRatio;
info.smoothness = smoothness;
end

function curve = average_obj_curve(objCell)
if isnumeric(objCell)
    curve = objCell(:);
    curve = curve(isfinite(curve));
    return;
end

if ~iscell(objCell)
    curve = [];
    return;
end

validCurves = {};
minLen = inf;
for i = 1:numel(objCell)
    obj = objCell{i};
    if isempty(obj)
        continue;
    end
    obj = obj(:);
    obj = obj(isfinite(obj));
    if isempty(obj)
        continue;
    end
    validCurves{end + 1} = obj; %#ok<AGROW>
    minLen = min(minLen, numel(obj));
end

if isempty(validCurves)
    curve = [];
    return;
end

M = zeros(minLen, numel(validCurves));
for i = 1:numel(validCurves)
    M(:, i) = validCurves{i}(1:minLen);
end

curve = mean(M, 2);
end

function [alphaVal, betaVal, lambdaVal, kVal] = parse_params_from_path(modelFile)
parts = split(string(modelFile), filesep);
alphaVal = NaN;
betaVal = NaN;
lambdaVal = NaN;
kVal = NaN;

for i = 1:numel(parts)
    item = parts(i);
    if startsWith(item, 'alpha-')
        alphaVal = str2double(extractAfter(item, 'alpha-'));
    elseif startsWith(item, 'beta-')
        betaVal = str2double(extractAfter(item, 'beta-'));
    elseif startsWith(item, 'lambda-')
        lambdaVal = str2double(extractAfter(item, 'lambda-'));
    elseif startsWith(item, 'k-')
        kVal = str2double(extractAfter(item, 'k-'));
    end
end
end
