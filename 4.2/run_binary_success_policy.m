function run_binary_success_policy(policyName, employedRule, employedOperator, onlookerRule, onlookerOperator)
scriptDir = fileparts(mfilename('fullpath'));
codeDir = fileparts(scriptDir);
addpath(fullfile(codeDir, 'common'));
fhd = str2func('cec14_func');
functions = [3, 5, 19];
NP = 50;
dimension = 30;
maxFEs = 10000 * dimension;
Limit = 50;
numRuns = 30;
c = 1;
W = 1;
for TF = functions
    zong_time = zeros(numRuns, 1);
    zong_zuiyou = zeros(maxFEs, numRuns);
    zong_success_history = zeros(maxFEs, numRuns);
    for run_i = 1:numRuns
        rng(420000000 + TF * 1000 + run_i, 'twister');
        objectiveFcn = @(x) feval(fhd, x', TF);
        cfg = struct('NP', NP, 'maxFEs', maxFEs, 'limit', Limit, ...
            'c', c, 'rewardWeight', W, 'beta', 1.5, ...
            'alphaLevy', 1, 'M', 1, 'initializationOrder', 'random', ...
            'employedRule', employedRule, 'onlookerRule', onlookerRule, ...
            'employedFixedOperator', employedOperator, ...
            'onlookerFixedOperator', onlookerOperator);
        timerValue = tic;
        runResult = mab_abc_core(objectiveFcn, dimension, -100, 100, cfg);
        zong_time(run_i) = toc(timerValue);
        zong_zuiyou(:, run_i) = runResult.bestHistory;
        zong_success_history(:, run_i) = runResult.cumulativeRewardHistory;
    end
    finalValues = zong_zuiyou(maxFEs, :);
    zong_zuiyou_pingjunzhi = mean(finalValues);
    zong_zuiyou_std = std(finalValues, 0);
    filename = sprintf('zong_ABC_MAB_FEs_%s_D%d_TF%d.mat', ...
        policyName, dimension, TF);
    save(filename, 'policyName', 'zong_time', 'zong_zuiyou', ...
        'zong_success_history', 'zong_zuiyou_pingjunzhi', ...
        'zong_zuiyou_std', 'NP', 'maxFEs', 'dimension', 'Limit', ...
        'TF', 'W', 'c', 'employedRule', 'employedOperator', ...
        'onlookerRule', 'onlookerOperator');
end
end
