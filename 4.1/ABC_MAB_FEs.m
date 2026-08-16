clc;
clearvars;
scriptDir = fileparts(mfilename('fullpath'));
codeDir = fileparts(scriptDir);
addpath(fullfile(codeDir, 'common'));
fhd = str2func('cec14_func');
NP = 50;
dimension = 30;
maxFEs = 10000 * dimension;
Limit = 50;
W = 1;
rhoValues = [0.1, 0.5, 1.0, 1.5, 2.0, 5.0];
numRuns = 50;
lb = -100;
ub = 100;
for rhoIndex = 1:numel(rhoValues)
    c = rhoValues(rhoIndex) * W;
    for TF = 1:30
        zong_time = zeros(numRuns, 1);
        zong_zuiyou = zeros(maxFEs, numRuns);
        for run_i = 1:numRuns
            rng(410000000 + rhoIndex * 100000 + TF * 1000 + run_i, 'twister');
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
        filename = sprintf('zong_ABC_MAB_FEs_D%d_TF%d_rho_%g.mat', ...
            dimension, TF, rhoValues(rhoIndex));
        save(filename, 'zong_time', 'zong_zuiyou', ...
            'zong_zuiyou_pingjunzhi', 'zong_zuiyou_std', 'NP', ...
            'maxFEs', 'dimension', 'Limit', 'TF', 'W', 'c');
    end
end
