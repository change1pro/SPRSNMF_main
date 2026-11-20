clear
clc

rootDir = fileparts(mfilename('fullpath'));
paramDir = fullfile(rootDir, 'SPRSNMF_top_params');
outDir = fullfile(rootDir, 'visualization_results', ...
    'alpha_lambda_sensitivity');
if exist(outDir, 'dir') == 0
    mkdir(outDir);
end

datasets = {'att','UMIST','JAFFE','oh15','GLIOMA', ...
    'BA','lymphoma','lung','colon'};
metricNames = {'ACC_mean','NMI_mean','Purity_mean'};
metricLabels = {'ACC','NMI','Purity'};
figureNames = {'SPRSNMF_alpha_lambda_ACC', ...
    'SPRSNMF_alpha_lambda_NMI', ...
    'SPRSNMF_alpha_lambda_Purity'};

% One color for each alpha value, from blue to yellow.
rowColors = [ ...
    0.30, 0.23, 0.71; ...
    0.20, 0.41, 0.91; ...
    0.14, 0.67, 0.96; ...
    0.14, 0.76, 0.73; ...
    0.34, 0.81, 0.48; ...
    0.96, 0.67, 0.23; ...
    0.96, 0.91, 0.18];

fixed_rows = cell(numel(datasets), 4);

for metricId = 1:numel(metricNames)
    fig = figure('Color', 'w', 'Position', [70, 55, 1280, 860]);
    columnLeft = [0.060, 0.355, 0.650];
    rowBottom = [0.655, 0.330, 0.005];
    axesWidth = 0.285;
    axesHeight = 0.340;
    for d = 1:numel(datasets)
        dataset = datasets{d};
        coarseFile = fullfile(paramDir, ...
            [dataset '_all_coarse_sorted.csv']);
        topFile = fullfile(paramDir, ...
            [dataset '_top50_params.csv']);

        if exist(coarseFile, 'file') ~= 2
            error('Missing coarse-search table: %s', coarseFile);
        end
        if exist(topFile, 'file') ~= 2
            error('Missing top-parameter table: %s', topFile);
        end

        T = readtable(coarseFile, 'VariableNamingRule', 'preserve');
        topT = readtable(topFile, 'VariableNamingRule', 'preserve');

        required = {'alpha','beta','lambda','k',metricNames{metricId}};
        if ~all(ismember(required, T.Properties.VariableNames))
            error('Required columns are missing from: %s', coarseFile);
        end
        if isempty(topT)
            error('Top-parameter table is empty: %s', topFile);
        end

        % Hold beta and k at the values of the best coarse-search setting.
        fixedBeta = topT.beta(1);
        fixedK = topT.k(1);
        betaTol = max(1e-12, abs(fixedBeta) * 1e-12);
        subT = T(abs(T.beta - fixedBeta) <= betaTol & T.k == fixedK, :);

        alphaVals = sort(unique(subT.alpha));
        lambdaVals = sort(unique(subT.lambda));
        if numel(alphaVals) ~= 7 || numel(lambdaVals) ~= 7
            error(['%s does not contain a complete 7-by-7 alpha-lambda ', ...
                'grid at beta=%g and k=%g.'], dataset, fixedBeta, fixedK);
        end

        Z = nan(numel(lambdaVals), numel(alphaVals));
        for ai = 1:numel(alphaVals)
            alphaTol = max(1e-12, abs(alphaVals(ai)) * 1e-12);
            for li = 1:numel(lambdaVals)
                lambdaTol = max(1e-12, abs(lambdaVals(li)) * 1e-12);
                row = subT(abs(subT.alpha - alphaVals(ai)) <= alphaTol & ...
                    abs(subT.lambda - lambdaVals(li)) <= lambdaTol, :);
                if height(row) ~= 1
                    error(['Expected one result for %s at alpha=%g, ', ...
                        'lambda=%g, beta=%g, k=%g; found %d.'], ...
                        dataset, alphaVals(ai), lambdaVals(li), ...
                        fixedBeta, fixedK, height(row));
                end
                Z(li, ai) = 100 * row.(metricNames{metricId})(1);
            end
        end

        if metricId == 1
            fixed_rows(d,:) = {dataset, fixedBeta, fixedK, height(subT)};
        end

        rowId = ceil(d / 3);
        columnId = mod(d - 1, 3) + 1;
        ax = axes(fig, ...
            'OuterPosition', [columnLeft(columnId), rowBottom(rowId), ...
            axesWidth, axesHeight], ...
            'PositionConstraint', 'outerposition');

        bars = bar3(ax, Z);
        for barId = 1:numel(bars)
            colorId = min(barId, size(rowColors, 1));
            bars(barId).FaceColor = rowColors(colorId,:);
            bars(barId).EdgeColor = [0.22, 0.22, 0.22];
            bars(barId).LineWidth = 0.28;
        end

        title(ax, sprintf('(%c) %s', char('a' + d - 1), dataset), ...
            'FontName', 'Times New Roman', 'FontSize', 12, ...
            'FontWeight', 'bold');
        xlabel(ax, '\alpha', 'FontName', 'Times New Roman', ...
            'FontSize', 11, 'FontWeight', 'bold');
        ylabel(ax, '\lambda', 'FontName', 'Times New Roman', ...
            'FontSize', 11, 'FontWeight', 'bold');
        zlabel(ax, metricLabels{metricId}, ...
            'FontName', 'Times New Roman', 'FontSize', 11, ...
            'FontWeight', 'bold');

        set(ax, 'XTick', 1:numel(alphaVals), ...
            'XTickLabel', make_power_labels(alphaVals), ...
            'YTick', 1:numel(lambdaVals), ...
            'YTickLabel', make_power_labels(lambdaVals), ...
            'FontName', 'Times New Roman', 'FontSize', 10, ...
            'FontWeight', 'bold', 'Color', 'w', 'LineWidth', 0.9);
        xtickangle(ax, 45);
        ytickangle(ax, 45);
        grid(ax, 'on');
        set(ax, 'GridColor', [0.84, 0.84, 0.84], ...
            'GridAlpha', 0.42);
        view(ax, -35, 28);
    end

    drawnow;
    baseName = fullfile(outDir, figureNames{metricId});
    exportgraphics(fig, [baseName '.png'], ...
        'Resolution', 300, 'BackgroundColor', 'white');
    exportgraphics(fig, [baseName '.pdf'], ...
        'ContentType', 'vector', 'BackgroundColor', 'white');
    savefig(fig, [baseName '.fig']);
    close(fig);
end

fixedTable = cell2table(fixed_rows, ...
    'VariableNames', {'Dataset','Fixed_beta','Fixed_k','Grid_rows'});
writetable(fixedTable, ...
    fullfile(outDir, 'SPRSNMF_alpha_lambda_fixed_parameters.csv'));

fprintf('Alpha-lambda sensitivity figures will be saved to:\n%s\n', outDir);

function labels = make_power_labels(values)
labels = cell(size(values));
for i = 1:numel(values)
    exponent = round(log10(values(i)));
    if abs(values(i) - 10^exponent) <= max(1e-12, abs(values(i))*1e-12)
        if exponent == 0
            labels{i} = '1';
        else
            labels{i} = sprintf('10^{%d}', exponent);
        end
    else
        labels{i} = num2str(values(i));
    end
    labels{i} = strrep(labels{i}, '-', char(8722));
end
end
