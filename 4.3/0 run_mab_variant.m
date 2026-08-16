function Zuiyou = run_mab_variant(variant_name, fhd, TF, NP, dimension, ...
    maxFEs, Limit, c, w, beta, alpha_levy, M, lb, ub)
    switch variant_name
        case 'Phase-UCB'
            employed_pool = [1, 2, 3];
            onlooker_pool = [1, 4];
            selection_rule = 'ucb';
            controller_structure = 'dual';
        case 'Phase-Random'
            employed_pool = [1, 2, 3];
            onlooker_pool = [1, 4];
            selection_rule = 'random';
            controller_structure = 'dual';
        case 'Global-UCB'
            employed_pool = 1:4;
            onlooker_pool = 1:4;
            selection_rule = 'ucb';
            controller_structure = 'global';
        case 'Global-Random'
            employed_pool = 1:4;
            onlooker_pool = 1:4;
            selection_rule = 'random';
            controller_structure = 'global';
        case 'Dual-All-UCB'
            employed_pool = 1:4;
            onlooker_pool = 1:4;
            selection_rule = 'ucb';
            controller_structure = 'dual';
        otherwise
            error('Unknown variant: %s', variant_name);
    end
    C = lb + rand(NP, dimension) .* (ub - lb);
    Fitness = inf(1, NP);
    trial = zeros(1, NP);
    FEs = 0;
    Optimalvalue = inf;
    Zuiyou = zeros(maxFEs, 1);
    for i = 1:NP
        Fitness(i) = feval(fhd, C(i, :)', TF);
        FEs = FEs + 1;
        if Fitness(i) < Optimalvalue
            Optimalvalue = Fitness(i);
        end
        Zuiyou(FEs) = Optimalvalue;
    end
    n_emp = zeros(1, numel(employed_pool));
    mean_emp = zeros(1, numel(employed_pool));
    t_emp = 0;
    n_onl = zeros(1, numel(onlooker_pool));
    mean_onl = zeros(1, numel(onlooker_pool));
    t_onl = 0;
    global_pool = 1:4;
    n_global = zeros(1, numel(global_pool));
    mean_global = zeros(1, numel(global_pool));
    t_global = 0;
    while FEs < maxFEs
        for i = 1:NP
            if FEs >= maxFEs
                break;
            end
            if strcmp(selection_rule, 'random')
                local_arm = randi(numel(employed_pool));
            elseif strcmp(controller_structure, 'global')
                t_global = t_global + 1;
                local_arm = select_ucb1_arm(n_global, mean_global, t_global, c);
            else
                t_emp = t_emp + 1;
                local_arm = select_ucb1_arm(n_emp, mean_emp, t_emp, c);
            end
            operator_id = employed_pool(local_arm);
            peer_index = random_peer(NP, i);
            C_new = apply_operator(C(i, :), C(peer_index, :), operator_id, ...
                dimension, beta, alpha_levy, FEs, M);
            C_new = max(min(C_new, ub), lb);
            fit_new = feval(fhd, C_new', TF);
            FEs = FEs + 1;
            reward = 0;
            if fit_new < Optimalvalue
                Optimalvalue = fit_new;
            end
            if fit_new < Fitness(i)
                C(i, :) = C_new;
                Fitness(i) = fit_new;
                trial(i) = 0;
                reward = w;
            else
                trial(i) = trial(i) + 1;
            end
            if strcmp(selection_rule, 'ucb')
                if strcmp(controller_structure, 'global')
                    n_global(local_arm) = n_global(local_arm) + 1;
                    mean_global(local_arm) = mean_global(local_arm) + ...
                        (reward - mean_global(local_arm)) / n_global(local_arm);
                else
                    n_emp(local_arm) = n_emp(local_arm) + 1;
                    mean_emp(local_arm) = mean_emp(local_arm) + ...
                        (reward - mean_emp(local_arm)) / n_emp(local_arm);
                end
            end
            Zuiyou(FEs) = Optimalvalue;
        end
        if FEs >= maxFEs
            break;
        end
        fitness_min = min(Fitness);
        quality = 1 ./ (Fitness - fitness_min + 1);
        probability = quality ./ sum(quality);
        for onlooker_i = 1:NP
            if FEs >= maxFEs
                break;
            end
            source_index = roulette_index(probability);
            peer_index = random_peer(NP, source_index);
            if strcmp(selection_rule, 'random')
                local_arm = randi(numel(onlooker_pool));
            elseif strcmp(controller_structure, 'global')
                t_global = t_global + 1;
                local_arm = select_ucb1_arm(n_global, mean_global, t_global, c);
            else
                t_onl = t_onl + 1;
                local_arm = select_ucb1_arm(n_onl, mean_onl, t_onl, c);
            end
            operator_id = onlooker_pool(local_arm);
            C_new = apply_operator(C(source_index, :), C(peer_index, :), ...
                operator_id, dimension, beta, alpha_levy, FEs, M);
            C_new = max(min(C_new, ub), lb);
            fit_new = feval(fhd, C_new', TF);
            FEs = FEs + 1;
            reward = 0;
            if fit_new < Optimalvalue
                Optimalvalue = fit_new;
            end
            if fit_new < Fitness(source_index)
                C(source_index, :) = C_new;
                Fitness(source_index) = fit_new;
                trial(source_index) = 0;
                reward = w;
            else
                trial(source_index) = trial(source_index) + 1;
            end
            if strcmp(selection_rule, 'ucb')
                if strcmp(controller_structure, 'global')
                    n_global(local_arm) = n_global(local_arm) + 1;
                    mean_global(local_arm) = mean_global(local_arm) + ...
                        (reward - mean_global(local_arm)) / n_global(local_arm);
                else
                    n_onl(local_arm) = n_onl(local_arm) + 1;
                    mean_onl(local_arm) = mean_onl(local_arm) + ...
                        (reward - mean_onl(local_arm)) / n_onl(local_arm);
                end
            end
            Zuiyou(FEs) = Optimalvalue;
        end
        if FEs >= maxFEs
            break;
        end
        for i = 1:NP
            if FEs >= maxFEs
                break;
            end
            if trial(i) > Limit
                C(i, :) = lb + rand(1, dimension) .* (ub - lb);
                Fitness(i) = feval(fhd, C(i, :)', TF);
                FEs = FEs + 1;
                trial(i) = 0;
                if Fitness(i) < Optimalvalue
                    Optimalvalue = Fitness(i);
                end
                Zuiyou(FEs) = Optimalvalue;
            end
        end
    end
    if FEs < maxFEs
        Zuiyou(FEs + 1:maxFEs) = Optimalvalue;
    end
end
function selected_arm = select_ucb1_arm(counts, mean_rewards, round_index, c)
    untried = find(counts == 0);
    if ~isempty(untried)
        selected_arm = untried(randi(numel(untried)));
        return;
    end
    ucb_scores = mean_rewards + c .* sqrt(log(round_index) ./ counts);
    best_score = max(ucb_scores);
    tolerance = 8 * eps(max(1, abs(best_score)));
    tied_arms = find(abs(ucb_scores - best_score) <= tolerance);
    selected_arm = tied_arms(randi(numel(tied_arms)));
end
function C_new = apply_operator(X_i, X_peer, operator_id, dim, beta, ...
    alpha_levy, FEs, M)
    C_new = X_i;
    switch operator_id
        case 1
            component = randi(dim);
            phi = 2 * rand - 1;
            C_new(component) = X_i(component) + ...
                phi * (X_i(component) - X_peer(component));
        case 2
            levy_peer = levy_perturb(X_peer, beta, alpha_levy);
            phi = 2 * rand(1, dim) - 1;
            C_new = X_i + phi .* (X_i - levy_peer);
        case 3
            levy_peer = levy_perturb(X_peer, beta, alpha_levy);
            phi = 2 * rand(1, dim) - 1;
            C_new = phi .* levy_peer;
        case 4
            component = randi(dim);
            phi = 2 * rand - 1;
            denominator = max(FEs * M, realmin);
            C_new(component) = X_i(component) + ...
                phi * (X_i(component) - X_peer(component)) / denominator;
        otherwise
            error('Unknown operator identifier.');
    end
end
function perturbed = levy_perturb(position, beta, alpha_levy)
    sigma_u = (gamma(1 + beta) * sin(pi * beta / 2) / ...
        (gamma((1 + beta) / 2) * beta * ...
        2^((beta - 1) / 2)))^(1 / beta);
    u = sigma_u .* randn(size(position));
    v = randn(size(position));
    levy_step = u ./ max(abs(v), realmin).^(1 / beta);
    perturbed = position + alpha_levy .* levy_step;
end
function peer_index = random_peer(population_size, excluded_index)
    peer_index = randi(population_size - 1);
    if peer_index >= excluded_index
        peer_index = peer_index + 1;
    end
end
function index = roulette_index(probability)
    threshold = rand;
    cumulative_probability = cumsum(probability);
    index = find(threshold <= cumulative_probability, 1, 'first');
    if isempty(index)
        index = numel(probability);
    end
end
