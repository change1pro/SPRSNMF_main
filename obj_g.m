clear
clc
addpath(genpath(pwd));
addpath(fullfile(pwd, 'utils'), '-begin');
addpath(fullfile(pwd, 'core'), '-begin');

data_names = {'Yale', 'COIL20', 'JAFFE', 'ORL', 'GLIOMA', 'BA', 'att', 'UMIST'};

alpha_list = [10^-6, 10^-4, 10^-2, 1, 10^2, 10^4, 10^6];
beta_list = [10^-6, 10^-4, 10^-2, 1, 10^2, 10^4, 10^6];
lambda_list = [10^-6, 10^-4, 10^-2, 1, 10^2, 10^4, 10^6];
k_list = [2, 3, 5, 7, 10, 15, 20];

obj_folder = 'obj';
if exist(obj_folder, 'dir') == 0
    mkdir(obj_folder);
end

for each_data = 1:numel(data_names)
    dataset = data_names{each_data};
    disp(dataset)

    [X, ~, info] = load_dataset(pwd, dataset);
    c = info.classes;

    data_folder = fullfile(obj_folder, dataset, 'all_params');
    if exist(data_folder, 'dir') == 0
        mkdir(data_folder);
    end

    for i = 1:numel(alpha_list)
        alpha = alpha_list(i);
        for j = 1:numel(beta_list)
            beta = beta_list(j);
            for t = 1:numel(lambda_list)
                lambda = lambda_list(t);
                for q = 1:numel(k_list)
                    k = min(k_list(q), info.samples - 1);

                    filename = fullfile(data_folder, ...
                        sprintf('alpha=%g_beta=%g_lambda=%g_k=%g.mat', alpha, beta, lambda, k));

                    if exist(filename, 'file') == 2
                        continue;
                    end

                    [~, obj] = splagnmf(X, c, alpha, lambda, beta, k);
                    save(filename, 'obj', 'alpha', 'beta', 'lambda', 'k');
                end
            end
        end
    end
end
