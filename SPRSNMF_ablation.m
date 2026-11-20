clear
clc

rootDir = fileparts(mfilename('fullpath'));
addpath(genpath(rootDir));
addpath(fullfile(rootDir, 'utils'), '-begin');
addpath(fullfile(rootDir, 'core'), '-begin');

quickTest = strcmp(getenv('SPLAGNMF_ABLATION_QUICK_TEST'), '1');

dataNames = {'att', 'UMIST', 'JAFFE', 'oh15', 'GLIOMA', ...
    'BA', 'lymphoma', 'lung', 'colon'};
parameterValues = [10^-6, 10^-4, 10^-2, 1, 10^2, 10^4, 10^6];
kValues = [3, 5, 7, 10, 15, 20];
coarseRuns = 3;
validationRuns = 10;
topNumber = 50;

if quickTest
    dataNames = {'GLIOMA'};
    parameterValues = 1;
    kValues = 3;
    coarseRuns = 1;
    validationRuns = 1;
    topNumber = 1;
    outputRoot = fullfile(tempdir, 'SPLAGNMF_ablation_search_quick_test');
else
    outputRoot = fullfile(rootDir, 'SPLAGNMF_ablation_search');
end

if exist(outputRoot, 'dir') == 0
    mkdir(outputRoot);
end

variants = struct( ...
    'key', {'without_graph', 'without_sparse', 'without_spl'}, ...
    'label', {'w/o Graph', 'w/o Sparse', 'w/o SPL'}, ...
    'method', {@splagnmf_without_graph, ...
               @splagnmf_without_sparse, ...
               @splagnmf_without_spl});

config.protocol = 'strict-independent-tuning-v1';
config.dataNames = dataNames;
config.parameterValues = parameterValues;
config.kValues = kValues;
config.coarseRuns = coarseRuns;
config.validationRuns = validationRuns;
config.topNumber = topNumber;
config.variantKeys = {variants.key};
config.quickTest = quickTest;

verify_or_save_config(outputRoot, config);

fprintf('Ablation output folder:\n%s\n', outputRoot);
fprintf('Quick-test mode: %d\n', quickTest);

for variantIndex = 1:numel(variants)
    variant = variants(variantIndex);
    parameterGrid = make_parameter_grid(variant.key, ...
        parameterValues, kValues);

    fprintf('\n============================================================\n');
    fprintf('Variant: %s | coarse combinations per dataset: %d\n', ...
        variant.label, height(parameterGrid));
    fprintf('============================================================\n');

    for datasetIndex = 1:numel(dataNames)
        dataset = dataNames{datasetIndex};
        fprintf('\nVariant %s | dataset %s\n', variant.label, dataset);

        [X, Y, info] = load_dataset(rootDir, dataset);
        c = info.classes;

        coarseFolder = fullfile(outputRoot, variant.key, ...
            'coarse', dataset);
        ensure_folder(coarseFolder);

        run_coarse_search(variant, parameterGrid, X, Y, info, c, ...
            datasetIndex, coarseRuns, coarseFolder);

        topFolder = fullfile(outputRoot, variant.key, 'top_params');
        ensure_folder(topFolder);
        topTable = select_top_parameters(variant, dataset, ...
            parameterGrid, coarseFolder, topFolder, topNumber);

        validationFolder = fullfile(outputRoot, variant.key, ...
            'validate', dataset);
        ensure_folder(validationFolder);
        run_validation(variant, topTable, X, Y, info, c, ...
            datasetIndex, validationRuns, validationFolder);

        write_combined_summary(rootDir, outputRoot, variants, dataNames);
    end
end

summaryFile = write_combined_summary(rootDir, outputRoot, ...
    variants, dataNames);

fprintf('\nAll requested ablation searches are complete.\n');
fprintf('Combined summary:\n%s\n', summaryFile);


function grid = make_parameter_grid(variantKey, parameterValues, kValues)
switch variantKey
    case 'without_graph'
        beta = parameterValues(:);
        rowCount = numel(beta);
        alpha = nan(rowCount, 1);
        lambda = nan(rowCount, 1);
        k = nan(rowCount, 1);

    case 'without_sparse'
        [alphaGrid, lambdaGrid, kGrid] = ndgrid( ...
            parameterValues, parameterValues, kValues);
        alpha = alphaGrid(:);
        beta = nan(numel(alpha), 1);
        lambda = lambdaGrid(:);
        k = kGrid(:);

    case 'without_spl'
        [alphaGrid, betaGrid, lambdaGrid, kGrid] = ndgrid( ...
            parameterValues, parameterValues, parameterValues, kValues);
        alpha = alphaGrid(:);
        beta = betaGrid(:);
        lambda = lambdaGrid(:);
        k = kGrid(:);

    otherwise
        error('Unknown ablation variant: %s', variantKey);
end

grid = table(alpha, beta, lambda, k);
end


function run_coarse_search(variant, parameterGrid, X, Y, info, c, ...
        datasetIndex, runs, coarseFolder)

combinationCount = height(parameterGrid);
for combinationIndex = 1:combinationCount
    settingFolder = fullfile(coarseFolder, ...
        sprintf('combo-%04d', combinationIndex));
    ensure_folder(settingFolder);

    modelFile = fullfile(settingFolder, 'model.mat');
    partialFile = fullfile(settingFolder, 'partial.mat');
    errorFile = fullfile(settingFolder, 'error.mat');

    if exist(modelFile, 'file') == 2
        continue;
    end

    params = parameterGrid(combinationIndex, :);
    fprintf('Coarse %s | combination %d/%d | %s\n', ...
        variant.label, combinationIndex, combinationCount, ...
        format_parameters(params));

    model = evaluate_setting(variant, params, X, Y, info, c, ...
        datasetIndex, combinationIndex, runs, 202609260, ...
        partialFile, errorFile, 'coarse');
    save(modelFile, 'model');

    if exist(partialFile, 'file') == 2
        delete(partialFile);
    end
end
end


function topTable = select_top_parameters(variant, dataset, ...
        parameterGrid, coarseFolder, topFolder, topNumber)

modelFiles = dir(fullfile(coarseFolder, 'combo-*', 'model.mat'));
if numel(modelFiles) ~= height(parameterGrid)
    error(['Coarse search for %s/%s is incomplete: expected %d models, ' ...
        'found %d.'], variant.key, dataset, height(parameterGrid), ...
        numel(modelFiles));
end

rows = cell(numel(modelFiles), 14);
validCount = 0;
for fileIndex = 1:numel(modelFiles)
    modelFile = fullfile(modelFiles(fileIndex).folder, ...
        modelFiles(fileIndex).name);
    saved = load(modelFile, 'model');
    if ~isfield(saved, 'model') || ~isfield(saved.model, 'opts')
        continue;
    end

    model = saved.model;
    validCount = validCount + 1;
    rows(validCount, :) = {dataset, model.opts.combinationIndex, ...
        model.opts.alpha, model.opts.beta, model.opts.lambda, model.opts.k, ...
        model.average_acc, model.std_acc, ...
        model.average_nmi, model.std_nmi, ...
        model.average_purity, model.std_purity, ...
        model.average_score, modelFile};
end

rows = rows(1:validCount, :);
if isempty(rows)
    error('No valid coarse-search models found for %s/%s.', ...
        variant.key, dataset);
end

summary = cell2table(rows, 'VariableNames', ...
    {'Dataset', 'CombinationIndex', 'alpha', 'beta', 'lambda', 'k', ...
     'ACC_mean', 'ACC_std', 'NMI_mean', 'NMI_std', ...
     'Purity_mean', 'Purity_std', 'Score_mean', 'SourceFile'});
summary = sortrows(summary, 'Score_mean', 'descend');

topCount = min(topNumber, height(summary));
topTable = summary(1:topCount, :);

allCsv = fullfile(topFolder, ...
    sprintf('%s_all_coarse_sorted.csv', dataset));
topCsv = fullfile(topFolder, ...
    sprintf('%s_top%d_params.csv', dataset, topCount));
topMat = fullfile(topFolder, ...
    sprintf('%s_top%d_params.mat', dataset, topCount));

writetable(summary, allCsv);
writetable(topTable, topCsv);
save(topMat, 'topTable');
end


function run_validation(variant, topTable, X, Y, info, c, ...
        datasetIndex, runs, validationFolder)

for rankIndex = 1:height(topTable)
    settingFolder = fullfile(validationFolder, ...
        sprintf('rank-%03d', rankIndex));
    ensure_folder(settingFolder);

    modelFile = fullfile(settingFolder, 'model.mat');
    partialFile = fullfile(settingFolder, 'partial.mat');
    errorFile = fullfile(settingFolder, 'error.mat');

    if exist(modelFile, 'file') == 2
        continue;
    end

    params = topTable(rankIndex, {'alpha', 'beta', 'lambda', 'k'});
    fprintf('Validate %s | rank %d/%d | %s\n', ...
        variant.label, rankIndex, height(topTable), ...
        format_parameters(params));

    model = evaluate_setting(variant, params, X, Y, info, c, ...
        datasetIndex, rankIndex, runs, 202609261, ...
        partialFile, errorFile, 'validation');
    model.opts.coarseRank = rankIndex;
    model.opts.coarseCombinationIndex = ...
        topTable.CombinationIndex(rankIndex);
    save(modelFile, 'model');

    if exist(partialFile, 'file') == 2
        delete(partialFile);
    end
end
end


function model = evaluate_setting(variant, params, X, Y, info, c, ...
        datasetIndex, itemIndex, runs, seedBase, partialFile, ...
        errorFile, stage)

tempAcc = nan(1, runs);
tempNmi = nan(1, runs);
tempPurity = nan(1, runs);
tempTime = nan(1, runs);
objCell = cell(1, runs);
doneRuns = false(1, runs);

if exist(partialFile, 'file') == 2
    saved = load(partialFile);
    required = {'tempAcc', 'tempNmi', 'tempPurity', ...
        'tempTime', 'objCell', 'doneRuns'};
    if all(isfield(saved, required)) && numel(saved.doneRuns) == runs
        tempAcc = saved.tempAcc;
        tempNmi = saved.tempNmi;
        tempPurity = saved.tempPurity;
        tempTime = saved.tempTime;
        objCell = saved.objCell;
        doneRuns = saved.doneRuns;
        fprintf('  Resume: %d/%d runs already complete.\n', ...
            sum(doneRuns), runs);
    end
end

[alphaCall, betaCall, lambdaCall, kCall] = ...
    callable_parameters(params, info.samples);
method = variant.method;

for runIndex = 1:runs
    if doneRuns(runIndex)
        continue;
    end

    seed = seedBase + datasetIndex * 10000000 + ...
        itemIndex * 100 + runIndex;
    rng(seed, 'twister');

    try
        runTimer = tic;
        [F, obj] = method(X, c, alphaCall, lambdaCall, ...
            betaCall, kCall);
        pred = cluster_embedding(F, c);
        metric = clustering_metrics(Y, pred);
        elapsed = toc(runTimer);

        if any(~isfinite(F(:))) || any(~isfinite(obj(:)))
            error('AblationSearch:NonfiniteOutput', ...
                'The algorithm returned NaN or Inf values.');
        end

        tempAcc(runIndex) = metric.ACC;
        tempNmi(runIndex) = metric.NMI;
        tempPurity(runIndex) = metric.Purity;
        tempTime(runIndex) = elapsed;
        objCell{runIndex} = obj;
        doneRuns(runIndex) = true;

        save(partialFile, 'tempAcc', 'tempNmi', 'tempPurity', ...
            'tempTime', 'objCell', 'doneRuns', 'params', ...
            'runs', 'stage');

        fprintf(['  run %02d/%02d: ACC %.4f, NMI %.4f, ' ...
            'Purity %.4f, time %.2fs\n'], runIndex, runs, ...
            tempAcc(runIndex), tempNmi(runIndex), ...
            tempPurity(runIndex), tempTime(runIndex));
    catch ME
        failure.identifier = ME.identifier;
        failure.message = ME.message;
        failure.stack = ME.stack;
        failure.variant = variant.key;
        failure.stage = stage;
        failure.datasetIndex = datasetIndex;
        failure.itemIndex = itemIndex;
        failure.runIndex = runIndex;
        failure.params = params;
        save(errorFile, 'failure');
        rethrow(ME);
    end
end

model.average_acc = mean(tempAcc, 'omitnan');
model.std_acc = std(tempAcc, 0, 'omitnan');
model.average_nmi = mean(tempNmi, 'omitnan');
model.std_nmi = std(tempNmi, 0, 'omitnan');
model.average_purity = mean(tempPurity, 'omitnan');
model.std_purity = std(tempPurity, 0, 'omitnan');
model.average_score = mean([model.average_acc, ...
    model.average_nmi, model.average_purity]);
model.raw_acc = tempAcc;
model.raw_nmi = tempNmi;
model.raw_purity = tempPurity;
model.raw_time = tempTime;
model.obj = objCell;
model.time = sum(tempTime, 'omitnan');
model.opts.variant = variant.key;
model.opts.stage = stage;
model.opts.alpha = params.alpha;
model.opts.beta = params.beta;
model.opts.lambda = params.lambda;
model.opts.k = params.k;
model.opts.runs = runs;
model.opts.combinationIndex = itemIndex;
end


function [alpha, beta, lambda, k] = callable_parameters(params, sampleCount)
alpha = params.alpha;
beta = params.beta;
lambda = params.lambda;
k = params.k;

if isnan(alpha)
    alpha = 0;
end
if isnan(beta)
    beta = 0;
end
if isnan(lambda)
    lambda = 0;
end
if isnan(k)
    k = 1;
end

k = min(max(round(k), 1), max(sampleCount - 1, 1));
end


function text = format_parameters(params)
text = sprintf('alpha=%s, beta=%s, lambda=%s, k=%s', ...
    format_value(params.alpha), format_value(params.beta), ...
    format_value(params.lambda), format_value(params.k));
end


function text = format_value(value)
if isnan(value)
    text = 'N/A';
else
    text = num2str(value, '%.6g');
end
end


function summaryFile = write_combined_summary(rootDir, outputRoot, ...
        variants, dataNames)

rows = collect_full_model_rows(rootDir, outputRoot, dataNames);

for variantIndex = 1:numel(variants)
    variant = variants(variantIndex);
    for datasetIndex = 1:numel(dataNames)
        dataset = dataNames{datasetIndex};
        validationFolder = fullfile(outputRoot, variant.key, ...
            'validate', dataset);
        modelFiles = dir(fullfile(validationFolder, 'rank-*', 'model.mat'));
        if isempty(modelFiles)
            continue;
        end

        bestModel = [];
        bestFile = '';
        for fileIndex = 1:numel(modelFiles)
            modelFile = fullfile(modelFiles(fileIndex).folder, ...
                modelFiles(fileIndex).name);
            saved = load(modelFile, 'model');
            if ~isfield(saved, 'model')
                continue;
            end
            if isempty(bestModel) || ...
                    saved.model.average_score > bestModel.average_score
                bestModel = saved.model;
                bestFile = modelFile;
            end
        end

        if isempty(bestModel)
            continue;
        end

        rows(end + 1, :) = {variant.label, dataset, ...
            bestModel.average_acc, bestModel.std_acc, ...
            bestModel.average_nmi, bestModel.std_nmi, ...
            bestModel.average_purity, bestModel.std_purity, ...
            bestModel.average_score, bestModel.opts.alpha, ...
            bestModel.opts.beta, bestModel.opts.lambda, ...
            bestModel.opts.k, bestFile}; %#ok<AGROW>
    end
end

summaryFile = fullfile(outputRoot, 'SPLAGNMF_ablation_best_summary.csv');
if isempty(rows)
    return;
end

summary = cell2table(rows, 'VariableNames', ...
    {'Variant', 'Dataset', 'ACC_mean', 'ACC_std', ...
     'NMI_mean', 'NMI_std', 'Purity_mean', 'Purity_std', ...
     'Score_mean', 'alpha', 'beta', 'lambda', 'k', 'ModelFile'});
summary = sortrows(summary, {'Dataset', 'Variant'});
writetable(summary, summaryFile);
end


function rows = collect_full_model_rows(rootDir, outputRoot, dataNames)
rows = {};
cacheFile = fullfile(outputRoot, 'Full_best_summary_cache.csv');
officialFile = fullfile(rootDir, 'SPLAGNMF_top_validate', ...
    'SPLAGNMF_top_validate_best.csv');

if exist(officialFile, 'file') == 2
    fullTable = readtable(officialFile, 'TextType', 'string');
    if ismember('BestBy', fullTable.Properties.VariableNames)
        fullTable = fullTable(strcmpi(string(fullTable.BestBy), 'Score'), :);
    end
    rows = full_table_to_rows(fullTable);
    return;
end

if exist(cacheFile, 'file') == 2
    fullTable = readtable(cacheFile, 'TextType', 'string');
    rows = full_table_to_rows(fullTable);
    return;
end

for datasetIndex = 1:numel(dataNames)
    dataset = dataNames{datasetIndex};
    validationFolder = fullfile(rootDir, 'SPLAGNMF_top_validate', dataset);
    modelFiles = dir(fullfile(validationFolder, '**', 'model.mat'));
    bestModel = [];
    bestFile = '';

    for fileIndex = 1:numel(modelFiles)
        modelFile = fullfile(modelFiles(fileIndex).folder, ...
            modelFiles(fileIndex).name);
        saved = load(modelFile, 'model');
        if ~isfield(saved, 'model')
            continue;
        end
        if isempty(bestModel) || ...
                saved.model.average_score > bestModel.average_score
            bestModel = saved.model;
            bestFile = modelFile;
        end
    end

    if isempty(bestModel)
        continue;
    end

    rows(end + 1, :) = {'Full', dataset, ...
        bestModel.average_acc, bestModel.std_acc, ...
        bestModel.average_nmi, bestModel.std_nmi, ...
        bestModel.average_purity, bestModel.std_purity, ...
        bestModel.average_score, bestModel.opts.alpha, ...
        bestModel.opts.beta, bestModel.opts.lambda, ...
        bestModel.opts.k, bestFile}; %#ok<AGROW>
end

if ~isempty(rows)
    fullTable = cell2table(rows, 'VariableNames', ...
        {'Variant', 'Dataset', 'ACC_mean', 'ACC_std', ...
         'NMI_mean', 'NMI_std', 'Purity_mean', 'Purity_std', ...
         'Score_mean', 'alpha', 'beta', 'lambda', 'k', 'ModelFile'});
    writetable(fullTable, cacheFile);
end
end


function rows = full_table_to_rows(fullTable)
rows = {};
for rowIndex = 1:height(fullTable)
    rows(end + 1, :) = {'Full', char(string(fullTable.Dataset(rowIndex))), ...
        fullTable.ACC_mean(rowIndex), fullTable.ACC_std(rowIndex), ...
        fullTable.NMI_mean(rowIndex), fullTable.NMI_std(rowIndex), ...
        fullTable.Purity_mean(rowIndex), fullTable.Purity_std(rowIndex), ...
        fullTable.Score_mean(rowIndex), fullTable.alpha(rowIndex), ...
        fullTable.beta(rowIndex), fullTable.lambda(rowIndex), ...
        fullTable.k(rowIndex), char(string(fullTable.ModelFile(rowIndex)))}; %#ok<AGROW>
end
end


function verify_or_save_config(outputRoot, config)
configFile = fullfile(outputRoot, 'experiment_config.mat');
if exist(configFile, 'file') == 2
    saved = load(configFile, 'config');
    if ~isfield(saved, 'config') || ~isequaln(saved.config, config)
        error(['The existing ablation result folder uses a different ' ...
            'configuration. Use a new output folder or restore the ' ...
            'original settings.']);
    end
else
    save(configFile, 'config');
end
end


function ensure_folder(folderPath)
if exist(folderPath, 'dir') == 0
    mkdir(folderPath);
end
end
