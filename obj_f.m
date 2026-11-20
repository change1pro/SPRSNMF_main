clear
clc
addpath(genpath(pwd));

data_names = {'att','UMIST','JAFFE','oh15','GLIOMA','BA', ...
    'lymphoma','lung','colon'};

obj_folder = 'obj';
img_folder = 'plot_objimg';
if exist(img_folder, 'dir') == 0
    mkdir(img_folder);
end

for each_data = 1:numel(data_names)
    dataset = data_names{each_data};

    data_folder = fullfile(obj_folder, dataset, 'all_params');
    if exist(data_folder, 'dir') == 0
        disp(['No obj folder: ', data_folder]);
        continue;
    end

    files = dir(fullfile(data_folder, '*.mat'));
    if isempty(files)
        disp(['No obj files: ', data_folder]);
        continue;
    end

    save_folder = fullfile(img_folder, dataset);
    if exist(save_folder, 'dir') == 0
        mkdir(save_folder);
    end

    for i = 1:numel(files)
        obj_filename = fullfile(files(i).folder, files(i).name);
        load(obj_filename, 'obj', 'alpha', 'beta', 'lambda', 'k');

        if ~exist('obj', 'var') || isempty(obj)
            continue;
        end

        obj = obj(:);
        obj = obj(isfinite(obj));
        if numel(obj) < 2
            continue;
        end

        param_name = sprintf('alpha=%g, beta=%g, lambda=%g, k=%g', alpha, beta, lambda, k);
        file_name = sprintf('alpha=%g_beta=%g_lambda=%g_k=%g.png', alpha, beta, lambda, k);

        fontsize = 20;
        figure1 = figure('Color', [1 1 1], 'Visible', 'off');
        axes('Parent', figure1, 'FontSize', fontsize, 'FontName', 'Times New Roman');

        x = 1:length(obj);
        xx = 1:0.01:length(obj);
        yy = interp1(x, obj, xx, 'cubic');

        plot(x, obj, '.', xx, yy, 'linewidth', 3);
        xlabel('Iteration number', 'FontSize', fontsize);
        ylabel('Objective function value', 'FontSize', fontsize);
        title({dataset; param_name}, 'Interpreter', 'none');
        set(gca, 'FontName', 'Times New Roman', 'FontSize', fontsize);
        grid on;

        saveas(gcf, fullfile(save_folder, file_name));
        close(gcf);

        clear obj alpha beta lambda k
    end
end
