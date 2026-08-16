function plot_convergence_all11()
%PLOT_CONVERGENCE_ALL11 Plot all 11 algorithms in the same four panels.
% Two versions are exported:
%   1) median curves only;
%   2) median curves with interquartile-range bands.

    root_dir = fileparts(mfilename('fullpath'));
    output_dir = fullfile(root_dir, 'generated_figures');
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

%     functions = [3, 11, 19, 26];  7 8 11 14 16£¬
    functions = [3, 14, 18, 26];
    function_titles = { ...
        '(a) F3: Unimodal', ...
        '(b) F14: Simple multimodal', ...
        '(c) F18: Hybrid', ...
        '(d) F26: Composition'};
    checkpoint_pct = [1, 5, 10, 20, 50, 75, 100];
    max_fes = 300000;
    checkpoint_fes = round(checkpoint_pct / 100 * max_fes);
    error_floor = 1e-12;
    data_file_pattern = 'TF%d_cols50.mat';
    data_variable_name = 'data';

    algorithms = struct( ...
        'name', { ...
            'MAB-ABC', 'jSO', 'MadDE', 'LSHADE-SPACMA', ...
            'MPA', 'EO', 'BWO', 'ABCLIS', 'RNSABC', ...
            'BABC', 'ABC-AS^3'}, ...
        'folder', { ...
            'zong_ABC_MAB_c10_D30', 'zong_jSO_FEs_D30', 'madde', ...
            'zong_LSHADE_SPACMA_FEs_D30', 'mpa', 'EO', 'bwo', ...
            'ABC_LIS', 'ABC_RNS', 'ABC_B', 'zong_ABC_AS3_FEs_D30'}, ...
        'color', { ...
            [0.82, 0.08, 0.08], ...  % MAB-ABC
            [0.00, 0.38, 0.72], ...  % jSO
            [0.95, 0.45, 0.05], ...  % MadDE
            [0.10, 0.58, 0.18], ...  % LSHADE-SPACMA
            [0.50, 0.15, 0.58], ...  % MPA
            [0.00, 0.66, 0.68], ...  % EO
            [0.28, 0.28, 0.28], ...  % BWO
            [0.22, 0.58, 0.90], ...  % ABCLIS
            [0.88, 0.18, 0.45], ...  % RNSABC
            [0.45, 0.68, 0.08], ...  % BABC
            [0.58, 0.36, 0.78]}, ... % ABC-AS3
        'line_style', { ...
            '-', '-', '--', '-.', ':', '--', '-.', ':', '--', '-.', ':'}, ...
        'marker', { ...
            'o', 's', '^', 'd', 'v', '>', '<', 'p', 'h', 'x', '+'});

    statistics = collect_statistics(root_dir, algorithms, functions, ...
        checkpoint_fes, error_floor, data_file_pattern, data_variable_name);

    draw_all11(output_dir, algorithms, statistics, functions, ...
        function_titles, checkpoint_pct, false, ...
        'convergence_all11_median_D30_cols50');

    draw_all11(output_dir, algorithms, statistics, functions, ...
        function_titles, checkpoint_pct, true, ...
        'convergence_all11_iqr_D30_cols50');

    fprintf('All-11 convergence figures saved in: %s\n', output_dir);
end

function statistics = collect_statistics(root_dir, algorithms, functions, ...
    checkpoint_fes, error_floor, data_file_pattern, data_variable_name)

    statistics = repmat(struct('median', [], 'q1', [], 'q3', [], 'n', 0), ...
        numel(functions), numel(algorithms));

    for function_index = 1:numel(functions)
        tf = functions(function_index);
        optimum = 100 * tf;
        for algorithm_index = 1:numel(algorithms)
            algorithm = algorithms(algorithm_index);
            input_file = fullfile(root_dir, algorithm.folder, ...
                sprintf(data_file_pattern, tf));
            assert(exist(input_file, 'file') == 2, ...
                'Missing convergence file: %s', input_file);

            loaded = load(input_file, data_variable_name);
            assert(isfield(loaded, data_variable_name), ...
                'Variable "%s" was not found in %s.', data_variable_name, input_file);
            trajectories = loaded.(data_variable_name);
            assert(size(trajectories, 1) >= max(checkpoint_fes), ...
                'Insufficient FE records in %s.', input_file);
            checkpoint_values = trajectories(checkpoint_fes, :);
            errors = max(checkpoint_values - optimum, error_floor);

            statistics(function_index, algorithm_index).median = ...
                median(errors, 2);
            statistics(function_index, algorithm_index).q1 = ...
                row_percentile(errors, 25);
            statistics(function_index, algorithm_index).q3 = ...
                row_percentile(errors, 75);
            statistics(function_index, algorithm_index).n = ...
                size(trajectories, 2);
        end
    end
end

function draw_all11(output_dir, algorithms, statistics, functions, ...
    function_titles, checkpoint_pct, show_iqr, output_name)

    figure_handle = figure('Color', 'w', 'Position', [60, 40, 1400, 940], ...
        'Name', output_name, 'NumberTitle', 'off');
    axes_positions = [ ...
        0.065, 0.585, 0.415, 0.340; ...
        0.555, 0.585, 0.415, 0.340; ...
        0.065, 0.205, 0.415, 0.340; ...
        0.555, 0.205, 0.415, 0.340];

    line_handles = gobjects(numel(algorithms), 1);
    legend_labels = cell(numel(algorithms), 1);

    for function_index = 1:numel(functions)
        axes_handle = axes('Parent', figure_handle, ...
            'Position', axes_positions(function_index, :)); %#ok<LAXES>
        hold(axes_handle, 'on');

        if show_iqr
            for algorithm_index = 1:numel(algorithms)
                algorithm = algorithms(algorithm_index);
                current_stats = statistics(function_index, algorithm_index);
                band_x = [checkpoint_pct, fliplr(checkpoint_pct)];
                band_y = [current_stats.q1', fliplr(current_stats.q3')];
                % Preblend with white so the EPS version remains light.
                band_color = 0.94 * [1, 1, 1] + 0.06 * algorithm.color;
                fill(axes_handle, band_x, band_y, band_color, ...
                    'EdgeColor', 'none', 'HandleVisibility', 'off');
            end
        end

        for algorithm_index = 1:numel(algorithms)
            algorithm = algorithms(algorithm_index);
            current_stats = statistics(function_index, algorithm_index);
            line_width = 1.25;
            marker_size = 4.0;
            if algorithm_index == 1
                line_width = 2.8;
                marker_size = 6.2;
            end
            current_line = plot(axes_handle, checkpoint_pct, ...
                current_stats.median, ...
                'Color', algorithm.color, ...
                'LineStyle', algorithm.line_style, ...
                'LineWidth', line_width, ...
                'Marker', algorithm.marker, ...
                'MarkerSize', marker_size, ...
                'MarkerFaceColor', 'w');

            if function_index == 1
                line_handles(algorithm_index) = current_line;
                legend_labels{algorithm_index} = sprintf('%s (n=%d)', ...
                    algorithm.name, current_stats.n);
            end
        end

        set(axes_handle, 'YScale', 'log', 'FontName', 'Times New Roman', ...
            'FontSize', 9.5, 'LineWidth', 0.8, 'TickDir', 'out', ...
            'XTick', checkpoint_pct, 'XLim', [0, 101]);
        grid(axes_handle, 'on');
        axes_handle.GridAlpha = 0.17;
        axes_handle.MinorGridAlpha = 0.07;
        if function_index > 2
            xlabel(axes_handle, 'Evaluation budget used (%)', ...
                'FontName', 'Times New Roman', 'FontSize', 10.5);
        end
        ylabel(axes_handle, 'Best-so-far error (log scale)', ...
            'FontName', 'Times New Roman', 'FontSize', 10.5);
        title(axes_handle, function_titles{function_index}, ...
            'FontName', 'Times New Roman', 'FontSize', 12, ...
            'FontWeight', 'bold');
        box(axes_handle, 'on');
    end

    legend_handle = legend(line_handles, legend_labels, ...
        'Orientation', 'horizontal', 'Location', 'none', ...
        'FontName', 'Times New Roman', 'FontSize', 8.5, ...
        'Interpreter', 'tex', 'Box', 'off');
    legend_handle.NumColumns = 6;
    legend_handle.Position = [0.045, 0.025, 0.91, 0.090];

    png_file = fullfile(output_dir, [output_name, '.png']);
    eps_file = fullfile(output_dir, [output_name, '.eps']);
    print(figure_handle, png_file, '-dpng', '-r300');
    print(figure_handle, eps_file, '-depsc', '-painters');
    close(figure_handle);
end

function percentile_values = row_percentile(values, percentile)
% Toolbox-independent linear-interpolation percentile along each row.
    sorted_values = sort(values, 2);
    sample_count = size(sorted_values, 2);
    position = 1 + (sample_count - 1) * percentile / 100;
    lower_index = floor(position);
    upper_index = ceil(position);
    weight = position - lower_index;
    percentile_values = (1 - weight) .* sorted_values(:, lower_index) + ...
        weight .* sorted_values(:, upper_index);
end
