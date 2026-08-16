clc;
clearvars;
scriptDir = fileparts(mfilename('fullpath'));
codeDir = fileparts(scriptDir);
addpath(fullfile(codeDir, 'common'));
fhd = str2func('cec14_func');
dimensions = [10, 30, 50];
NP = 50;
Limit = 50;
numRuns = 50;
c = 1;
W = 1;
lb = -100;
ub = 100;
for dimension = dimensions
    maxFEs = 10000 * dimension;
    for TF = 1:30
        zong_time = zeros(numRuns, 1);
        zong_zuiyou = zeros(maxFEs, numRuns);
        for run_i = 1:numRuns
            rng(450000000 + dimension * 100000 + TF * 1000 + run_i, 'twister');
            objectiveFcn = @(x) feval(fhd, x', TF);
            cfg = struct('NP', NP, 'maxFEs', maxFEs, 'limit', Limit, ...
                'c', c, 'rewardWeight', W, 'beta', 1.5, ...
                'alphaLevy', 1, 'M', 1, 'initializationOrder', 'random');
            timerValue = tic;
            runResult = mab_abc_core(objectiveFcn, dimension, lb, ub, cfg);
            zong_time(run_i) = toc(timerValue);
            zong_zuiyou(:, run_i) = runResult.bestHistory;
        end
        finalValues = zong_zuiyou(maxFEs, :);
        zong_zuiyou_pingjunzhi = mean(finalValues);
        zong_zuiyou_std = std(finalValues, 0);
        filename = sprintf('zong_ABC_MAB_c10_D%d_TF%d.mat', dimension, TF);
        save(filename, 'zong_time', 'zong_zuiyou', ...
            'zong_zuiyou_pingjunzhi', 'zong_zuiyou_std', 'NP', ...
            'maxFEs', 'Limit', 'TF', 'dimension', 'c', 'W');
    end
end
