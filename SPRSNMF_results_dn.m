clear
clc

rootDir = fileparts(mfilename('fullpath'));
addpath(genpath(rootDir));
addpath(fullfile(rootDir, 'utils'), '-begin');
addpath(fullfile(rootDir, 'core'), '-begin');

% data_names = {'att','UMIST','JAFFE','oh15','GLIOMA','BA','lymphoma','lung','colon','COIL20','warpPIE10P','TOX_171'};
data_names = {'YaleB_32x32'};
%'TOX_171'data_names = {'TOX_171'};

 alpha_list = [10^-6, 10^-4, 10^-2, 1, 10^2, 10^4, 10^6];
 beta_list = [10^-6, 10^-4, 10^-2, 1, 10^2, 10^4, 10^6];
 lambda_list = [10^-6, 10^-4, 10^-2, 1, 10^2, 10^4, 10^6];
 k_list = [3, 5, 7, 10, 15, 20];
 runs = 3;




result_folder = fullfile(rootDir, 'SPLAGNMF_results_coarse');
if exist(result_folder, 'dir') == 0
    mkdir(result_folder);
end

combo = [];
for i = 1:numel(alpha_list)
    for j = 1:numel(beta_list)
        for t = 1:numel(lambda_list)
            for q = 1:numel(k_list)
                combo = [combo; i, j, t, q]; %#ok<AGROW>
            end
        end
    end
end
num_combo = size(combo, 1);

for each_data = 1:numel(data_names)
    dataset = data_names{each_data};
    fprintf('\nDataset: %s\n', dataset);

    dataset_folder = fullfile(result_folder, dataset);
    if exist(dataset_folder, 'dir') == 0
        mkdir(dataset_folder);
    end

    progress_file = fullfile(dataset_folder, 'coarse_progress.mat');
    start_combo = 1;

    if exist(progress_file, 'file') == 2
        s = load(progress_file, 'progress');
        if isfield(s, 'progress') && ...
                isfield(s.progress, 'num_combo') && ...
                s.progress.num_combo == num_combo
            start_combo = s.progress.combo_id;
            fprintf('Resume %s coarse search from combo %d/%d\n', ...
                dataset, start_combo, num_combo);
        else
            fprintf('Old progress format found for %s; start from combo 1.\n', dataset);
        end
    end

    if start_combo > num_combo
        fprintf('%s coarse search is already complete.\n', dataset);
        continue;
    end

    [X, Y, info] = load_dataset(rootDir, dataset);
    c = info.classes;

    for combo_id = start_combo:num_combo
        i = combo(combo_id, 1);
        j = combo(combo_id, 2);
        t = combo(combo_id, 3);
        q = combo(combo_id, 4);

        alpha = alpha_list(i);
        beta = beta_list(j);
        lambda = lambda_list(t);
        k = min(k_list(q), info.samples - 1);

        outPath = fullfile(dataset_folder, ...
            ['alpha-' num2str(alpha)], ...
            ['beta-' num2str(beta)], ...
            ['lambda-' num2str(lambda)], ...
            ['k-' num2str(k)]);

        if exist(outPath, 'dir') == 0
            mkdir(outPath);
        end

        file_name = fullfile(outPath, 'model.mat');

        if exist(file_name, 'file') == 2
            progress.combo_id = combo_id + 1;
            progress.num_combo = num_combo;
            progress.dataset = dataset;
            progress.status = 'next';
            save(progress_file, 'progress');
            continue;
        end

        progress.combo_id = combo_id;
        progress.num_combo = num_combo;
        progress.dataset = dataset;
        progress.alpha = alpha;
        progress.beta = beta;
        progress.lambda = lambda;
        progress.k = k;
        progress.status = 'running';
        save(progress_file, 'progress');

        fprintf('Coarse: %s | combo %d/%d | alpha=%g, beta=%g, lambda=%g, k=%g\n', ...
            dataset, combo_id, num_combo, alpha, beta, lambda, k);

        temp_acc = zeros(1, runs);
        temp_nmi = zeros(1, runs);
        temp_purity = zeros(1, runs);
        obj_cell = cell(1, runs);

        tic;
        for r = 1:runs
            rng(20260610 + each_data * 100000 + ...
                i * 10000 + j * 1000 + t * 100 + q * 10 + r);

            [F, obj] = splagnmf(X, c, alpha, lambda, beta, k);
            pred = cluster_embedding(F, c);
            metric = clustering_metrics(Y, pred);

            temp_acc(r) = metric.ACC;
            temp_nmi(r) = metric.NMI;
            temp_purity(r) = metric.Purity;
            obj_cell{r} = obj;

            fprintf('  run %02d/%02d: ACC %.4f, NMI %.4f, Purity %.4f\n', ...
                r, runs, temp_acc(r), temp_nmi(r), temp_purity(r));
        end
        time = toc;

        model.average_acc = mean(temp_acc);
        model.std_acc = std(temp_acc);
        model.average_nmi = mean(temp_nmi);
        model.std_nmi = std(temp_nmi);
        model.average_purity = mean(temp_purity);
        model.std_purity = std(temp_purity);
        model.average_score = mean([model.average_acc, ...
            model.average_nmi, model.average_purity]);
        model.raw_acc = temp_acc;
        model.raw_nmi = temp_nmi;
        model.raw_purity = temp_purity;
        model.obj = obj_cell;
        model.time = time;
        model.opts.alpha = alpha;
        model.opts.beta = beta;
        model.opts.lambda = lambda;
        model.opts.k = k;
        model.opts.runs = runs;

        save(file_name, 'model');

        progress.combo_id = combo_id + 1;
        progress.num_combo = num_combo;
        progress.dataset = dataset;
        progress.status = 'next';
        save(progress_file, 'progress');
    end

    progress.combo_id = num_combo + 1;
    progress.num_combo = num_combo;
    progress.dataset = dataset;
    progress.status = 'complete';
    save(progress_file, 'progress');

    fprintf('%s coarse search finished.\n', dataset);
end

fprintf('\nSPLAGNMF coarse search finished for all selected datasets.\n');
