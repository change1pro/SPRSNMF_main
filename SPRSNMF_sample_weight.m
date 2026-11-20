clear
clc

rootDir = fileparts(mfilename('fullpath'));
outDir = fullfile(rootDir, 'visualization_results', ...
    'selfpaced_growth_ratio_yang_fig9');
summaryFile = fullfile(outDir, ...
    'SPRSNMF_yang_fig9_growth_ratio_summary.xlsx');

if exist(summaryFile, 'file') ~= 2
    error('Summary table not found: %s', summaryFile);
end


summary = readtable(summaryFile, 'VariableNamingRule', 'preserve');
datasets = {'att', 'UMIST', 'JAFFE', 'oh15', 'GLIOMA', ...
    'BA', 'lymphoma', 'lung', 'colon'};
percentSteps = [10, 15, 25, 40];

requiredVariables = {'Dataset', 'PercentStep', ...
    'ACC_mean', 'NMI_mean', 'Purity_mean'};
missingVariables = setdiff(requiredVariables, ...
    summary.Properties.VariableNames);
if ~isempty(missingVariables)
    error('Missing columns in summary table: %s', ...
        strjoin(missingVariables, ', '));
end


numericVariables = {'PercentStep', 'ACC_mean', 'NMI_mean', 'Purity_mean'};
for i = 1:numel(numericVariables)
    variableName = numericVariables{i};
    summary.(variableName) = to_numeric_column( ...
        summary.(variableName), variableName);
end

for d = 1:numel(datasets)
    for p = 1:numel(percentSteps)
        rowExists = any(strcmpi(strtrim(string(summary.Dataset)), ...
            datasets{d}) & summary.PercentStep == percentSteps(p));
        if ~rowExists
            error('Missing summary row: dataset=%s, PercentStep=%g.', ...
                datasets{d}, percentSteps(p));
        end
    end
end

plot_all_datasets(summary, datasets, percentSteps, 'ACC', ...
    fullfile(outDir, 'SPRSNMF_yang_fig9_all_datasets_ACC_summary_only'));
plot_all_datasets(summary, datasets, percentSteps, 'NMI', ...
    fullfile(outDir, 'SPRSNMF_yang_fig9_all_datasets_NMI_summary_only'));
plot_all_datasets(summary, datasets, percentSteps, 'Purity', ...
    fullfile(outDir, 'SPRSNMF_yang_fig9_all_datasets_Purity_summary_only'));

fprintf('Figures regenerated directly from:\n%s\n', summaryFile);


function plot_all_datasets(summary, datasets, percentSteps, ...
    metricName, outBase)
fig = figure('Color', 'w', 'Position', [55, 35, 1280, 900]);
layout = tiledlayout(3, 3, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

for d = 1:numel(datasets)
    datasetRows = summary(strcmpi(strtrim(string(summary.Dataset)), ...
        datasets{d}), :);
    y = nan(size(percentSteps));

    for p = 1:numel(percentSteps)
        row = datasetRows(datasetRows.PercentStep == percentSteps(p), :);
        values = row.(meanField);
        y(p) = values(1);
    end

    ax = nexttile(layout);
    plot(ax, 1:numel(percentSteps), y, '-o', ...
        'Color', [1.00, 0.20, 0.16], 'LineWidth', 1.35, ...
        'MarkerSize', 5.0, 'MarkerFaceColor', 'w');
    title(ax, sprintf('(%c) %s', char(96 + d), datasets{d}), ...
        'FontName', 'Times New Roman', 'FontWeight', 'bold', ...
        'FontSize', 10.5);
    xlabel(ax, 'The setting ratio', 'FontName', 'Times New Roman', ...
        'FontWeight', 'bold', 'FontSize', 9.5);
    ylabel(ax, metricName, 'FontName', 'Times New Roman', ...
        'FontWeight', 'bold', 'FontSize', 9.5);
    set(ax, 'XTick', 1:numel(percentSteps), ...
        'XTickLabel', strcat(string(percentSteps), '%'), ...
        'FontName', 'Times New Roman', 'FontSize', 8.5, ...
        'LineWidth', 0.8);
    xlim(ax, [0.8, numel(percentSteps) + 0.2]);
    set_tight_limits(ax, y);
    format_y_ticks_at_most_two_decimals(ax);
    grid(ax, 'on');
    ax.GridAlpha = 0.12;
    box(ax, 'on');
    ax.Toolbar.Visible = 'off';
    disableDefaultInteractivity(ax);
end

save_figure(fig, outBase);
end


function set_tight_limits(ax, values)
values = values(isfinite(values));
if isempty(values)
    ylim(ax, [0, 1]);
    return;
end

valueMin = min(values);
valueMax = max(values);
if abs(valueMax - valueMin) < 1e-10
    padding = max(0.01, 0.04 * max(valueMax, 0.1));
else
    padding = max(0.005, 0.15 * (valueMax - valueMin));
end

lowerLimit = max(0, valueMin - padding);
upperLimit = min(1, valueMax + padding);
if upperLimit <= lowerLimit
    upperLimit = min(1, lowerLimit + 0.05);
end
ylim(ax, [lowerLimit, upperLimit]);
end


function format_y_ticks_at_most_two_decimals(ax)
rawLimits = ax.YLim;
rawRange = rawLimits(2) - rawLimits(1);
tickStep = max(0.01, ceil(rawRange / 7 * 100) / 100);
lowerLimit = floor(rawLimits(1) / tickStep) * tickStep;
upperLimit = ceil(rawLimits(2) / tickStep) * tickStep;
lowerLimit = max(0, lowerLimit);
upperLimit = min(1, upperLimit);
if upperLimit <= lowerLimit
    upperLimit = min(1, lowerLimit + tickStep);
end
tickValues = lowerLimit:tickStep:upperLimit;
ax.YLim = [lowerLimit, upperLimit];
ax.YTick = tickValues;
tickLabels = strings(size(tickValues));
for i = 1:numel(tickValues)
    label = sprintf('%.2f', tickValues(i));
    label = regexprep(label, '0+$', '');
    label = regexprep(label, '\.$', '');
    tickLabels(i) = string(label);
end
ax.YTickLabel = tickLabels;
end


function values = to_numeric_column(values, variableName)
if isnumeric(values)
    values = double(values);
else
    values = str2double(strtrim(string(values)));
end
if any(~isfinite(values))
    error('Column %s contains a nonnumeric or missing value.', variableName);
end
values = values(:);
end


function save_figure(fig, outBase)
savefig(fig, [outBase '.fig']);
exportgraphics(fig, [outBase '.png'], 'Resolution', 300);
try
    exportgraphics(fig, [outBase '.pdf'], 'ContentType', 'vector');
catch
    print(fig, [outBase '.pdf'], '-dpdf', '-bestfit');
end
close(fig);
end
