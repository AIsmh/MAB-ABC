function run_control_experiment(variantName, fileTag)
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
fhd = str2func('cec14_func');
NP = 50;
dimension = 30;
maxFEs = 10000 * dimension;
Limit = 50;
numRuns = 50;
c = 1;
W = 1;
for TF = 1:30
    zong_time = zeros(numRuns, 1);
    zong_zuiyou = zeros(maxFEs, numRuns);
    for run_i = 1:numRuns
        rng(430000000 + TF * 1000 + run_i, 'twister');
        timerValue = tic;
        zong_zuiyou(:, run_i) = run_mab_variant(variantName, fhd, TF, ...
            NP, dimension, maxFEs, Limit, c, W, 1.5, 1, 1, -100, 100);
        zong_time(run_i) = toc(timerValue);
    end
    finalValues = zong_zuiyou(maxFEs, :);
    zong_zuiyou_pingjunzhi = mean(finalValues);
    zong_zuiyou_std = std(finalValues, 0);
    filename = sprintf('zong_ABC_MAB_%s_D%d_TF%d.mat', fileTag, dimension, TF);
    save(filename, 'variantName', 'zong_time', 'zong_zuiyou', ...
        'zong_zuiyou_pingjunzhi', 'zong_zuiyou_std', 'NP', ...
        'maxFEs', 'Limit', 'TF', 'dimension', 'c', 'W');
end
end
