clear
clc
close all

rootDir = fileparts(mfilename('fullpath'));
result_folder = fullfile(rootDir, 'visualization_results', ...
    'robustness_gaussian_variance_rnmf_splagnmf_v2');
summary_file = fullfile(result_folder, ...
    'RNMF_SPLAGNMF_gaussian_summary.xlsx');

data_names = {'att', 'UMIST', 'JAFFE'};
method_names = {'RNMF', 'SPLAGNMF'};
display_names = {'RNMF', 'SPRSNMF'};
metric_columns = {'ACCMean', 'NMIMean', 'PurityMean'};
metric_names = {'ACC', 'NMI', 'Purity'};

if exist(summary_file, 'file') ~= 2
    error('Summary table not found: %s', summary_file);
end

summary_table = readtable(summary_file, ...
    'VariableNamingRule', 'preserve');
required_columns = {'Dataset', 'Method', 'NoiseVariance', ...
    metric_columns{:}};
missing_columns = setdiff(required_columns, ...
    summary_table.Properties.VariableNames);
if ~isempty(missing_columns)
    error('Summary table is missing required columns: %s', ...
        strjoin(missing_columns, ', '));
end

noise_column = str2double(string(summary_table.NoiseVariance));
noise_variances = sort(unique(noise_column(isfinite(noise_column))))';
plot_noise_variances = noise_variances(noise_variances > 0);
if isempty(plot_noise_variances)
    error('No positive noise variance was found in the summary table.');
end

metric_data = read_plot_data(summary_table, data_names, ...
    method_names, plot_noise_variances, metric_columns);

colors = [0.12, 0.43, 0.70; 0.82, 0.25, 0.18];
markers = {'d', 'o'};
plot_x = 1:numel(plot_noise_variances);
n_data = numel(data_names);
n_methods = numel(method_names);

fig = figure('Color', 'w', 'Position', [60, 40, 1280, 960]);
tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
legend_handles = gobjects(1, n_methods);
panel_index = 0;

for data_index = 1:n_data
    for metric_index = 1:numel(metric_names)
        values = metric_data{metric_index};
        panel_index = panel_index + 1;
        ax = nexttile;
        hold(ax, 'on');

        for method_index = 1:n_methods
            y = 100 * squeeze(values(data_index, method_index, :));
            curve = plot(ax, plot_x, y, '-', ...
                'Color', colors(method_index, :), ...
                'LineWidth', 1.8, ...
                'Marker', markers{method_index}, ...
                'MarkerSize', 5.5, ...
                'MarkerFaceColor', colors(method_index, :));
            if panel_index == 1
                legend_handles(method_index) = curve;
            end
        end

        title(ax, sprintf('(%c) %s', ...
            char(96 + panel_index), data_names{data_index}), ...
            'FontName', 'Times New Roman', ...
            'FontSize', 12, 'FontWeight', 'bold');
        xlabel(ax, 'Noise variance', ...
            'FontName', 'Times New Roman', ...
            'FontSize', 11, 'FontWeight', 'bold');
        ylabel(ax, [metric_names{metric_index} ' (%)'], ...
            'FontName', 'Times New Roman', ...
            'FontSize', 11, 'FontWeight', 'bold');
        xlim(ax, [plot_x(1), plot_x(end)]);
        xticks(ax, plot_x);
        xticklabels(ax, compose('%g', plot_noise_variances));
        grid(ax, 'on');
        box(ax, 'on');
        set(ax, 'FontName', 'Times New Roman', ...
            'FontSize', 10, 'FontWeight', 'bold', ...
            'LineWidth', 0.9, ...
            'GridColor', [0.84, 0.84, 0.84], ...
            'GridAlpha', 0.30);
        ax.Toolbar.Visible = 'off';
        disableDefaultInteractivity(ax);
    end
end

lgd = legend(legend_handles, display_names, ...
    'Orientation', 'horizontal', 'Box', 'off');
lgd.Layout.Tile = 'north';
lgd.FontName = 'Times New Roman';
lgd.FontSize = 11;

figure_name = 'RNMF_SPRSNMF_gaussian_3datasets_3metrics';
png_file = fullfile(result_folder, [figure_name '.png']);
pdf_file = fullfile(result_folder, [figure_name '.pdf']);
fig_file = fullfile(result_folder, [figure_name '.fig']);

exportgraphics(fig, png_file, 'Resolution', 300);
try
    exportgraphics(fig, pdf_file, 'ContentType', 'vector');
catch
    print(fig, pdf_file, '-dpdf', '-bestfit');
end
savefig(fig, fig_file);

fprintf('\nFinished. The 3-by-3 figure was saved to:\n%s\n', ...
    result_folder);

function metric_data = read_plot_data(summary_table, data_names, ...
        method_names, noise_variances, metric_columns)
n_data = numel(data_names);
n_methods = numel(method_names);
n_noise = numel(noise_variances);
n_metrics = numel(metric_columns);

metric_data = cell(1, n_metrics);
for metric_index = 1:n_metrics
    metric_data{metric_index} = nan(n_data, n_methods, n_noise);
end

dataset_column = strtrim(string(summary_table.Dataset));
method_column = strtrim(string(summary_table.Method));
noise_column = str2double(string(summary_table.NoiseVariance));

for data_index = 1:n_data
    dataset_match = strcmpi(dataset_column, data_names{data_index});
    for method_index = 1:n_methods
        if strcmpi(method_names{method_index}, 'SPLAGNMF')
            method_match = strcmpi(method_column, 'SPLAGNMF') ...
                | strcmpi(method_column, 'SPRSNMF');
        else
            method_match = strcmpi(method_column, ...
                method_names{method_index});
        end

        for noise_index = 1:n_noise
            tolerance = max(1e-12, ...
                10 * eps(max(1, abs(noise_variances(noise_index)))));
            noise_match = abs(noise_column ...
                - noise_variances(noise_index)) <= tolerance;
            row = find(dataset_match & method_match & noise_match);

            if numel(row) ~= 1
                error(['Expected exactly one row for dataset %s, ' ...
                    'method %s, noise variance %g, but found %d.'], ...
                    data_names{data_index}, method_names{method_index}, ...
                    noise_variances(noise_index), numel(row));
            end

            for metric_index = 1:n_metrics
                column = str2double(string( ...
                    summary_table.(metric_columns{metric_index})));
                value = column(row);
                if ~isfinite(value)
                    error(['Non-finite %s value for dataset %s, ' ...
                        'method %s, noise variance %g.'], ...
                        metric_columns{metric_index}, ...
                        data_names{data_index}, ...
                        method_names{method_index}, ...
                        noise_variances(noise_index));
                end
                metric_data{metric_index}(data_index, ...
                    method_index, noise_index) = value;
            end
        end
    end
end
end
