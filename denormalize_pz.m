%{
dataset = 'laubLoomis';
start_idx = 1;
end_idx = 2000;
% Squared Laub stuff========================
%{
dirs = {
    [0, 0, 0, 0, 0, 0, 1]; 
    [0, 0, 0, 0, 1, 0, 0];
    [1, 0, 0, 0, 0, 0, 0];
    [0, -1, 0, 0, 0, 0, 0];
    [0, 0, -1, 0, 0, 0, 0];
    [0, 0, 0, 0, 0, 1, 0];
};
bs = {
    0.137; 
    0.0685; 
    0.457;
    -1.354;
    -1.6505;
    0.045;
};
power_case = 2;
%}
% 31 Cases ============================
%{
dirs = {
    [0, 0, 0, 0, 0, 0, 1]; 
    [0, 0, 0, 0, 1, 0, 0];
    [1, 0, 0, 0, 0, 0, 0];
    [0, -1, 0, 0, 0, 0, 0];
    [0, 0, -1, 0, 0, 0, 0];
    [0, 0, 0, 0, 0, 1, 0];
};
bs = {
    0.105;
    0.06;
    0.44;
    -1.373;
    -1.6508;
    -0.0501;
};
power_case = 31;
%}
%}
%===========================================
%{
dataset = 'VanDelPol';
start_idx = 1;
end_idx = 1348;
% 31 degree raised van ==========================
%{
dirs = {
    [1, 0];
    [-1, 0];
    [0, 1];
    [0, -1];
};

bs = {
    -2.02594;
    -2.13816;
    -2.7183;
    -2.814;
};
power_case = 31;
%}
% Experiment for squared van ==================================
%{
dirs = {
    [1, 0];
    [-1, 0];
    [0, 1];
    [0, -1];
};

bs = {
    -2.0165;
    -2.138;
    -2.73;
    -2.804;
};
power_case = 2;
%}
%}

% Some general parameters setup
exp_num = length(dirs);
split_order = 40;

% first store time, second store memory, thrid store init_density, forth
% store max_density
result_mat = nan(end_idx, exp_num, 4); 

for i = start_idx:end_idx
    file_path_E  = fullfile(dataset, sprintf('E_interval_%d.mat', i));
    file_path_G  = fullfile(dataset, sprintf('G_interval_%d.mat', i));
    file_path_c  = fullfile(dataset, sprintf('c_interval_%d.mat', i));
    file_path_GI = fullfile(dataset, sprintf('GI_interval_%d.mat', i));
    [G, E, center, GI, alpha_size] = preprocess_data(file_path_E, file_path_G, file_path_c, file_path_GI, power_case);
    
    for j = 1:exp_num
        % Make up the direction and initialzied info
        exp_info = {dirs{j}, bs{j}, split_order};
        [a, b, adjusted_vector, adjusted_value, c, b_val, max_depth] = prepare_problem(exp_info, alpha_size, GI, center);
        
        % Record the memeory and density of the pZ
        init_mem = numel(G) + numel(E) + numel(GI);
        mem_track = MemTracker(init_mem);

        [temp_density, ~] = size(E);
        density_info = DensityInfo(temp_density);
        fprintf("Init density info is: %d\n", density_info.get_init());
        
        % Time the intersection checking.
        tStart = tic;
        check_result = depz_intersection(G, a, b, E, c, b_val, adjusted_vector, adjusted_value, ...
                                         mem_track, 0, max_depth, density_info);
        time_usage = toc(tStart);
        print_result(check_result, i, j, time_usage, mem_track);
        result_mat(i, j, 1) = time_usage;
        result_mat(i, j, 2) = mem_track.get();
        result_mat(i, j, 3) = density_info.get_init();
        result_mat(i, j, 4) = density_info.get_max();
    end
end
% Note print_result will panic if there's inconclusive result or intersect
total_time_per_exp = sum(result_mat(:, :, 1), 1, 'omitnan');
max_memory_per_exp = max(result_mat(:, :, 2), [], 1, 'omitnan');
%{
min_init_density_over_sets = min(result_mat(:, :, 3), [], 1, 'omitnan');
max_init_density_over_sets = max(result_mat(:, :, 3), [], 1, 'omitnan');
mean_init_density_over_sets = mean(result_mat(:, :, 3), 1, 'omitnan');
min_max_density_over_sets = min(result_mat(:, :, 4), [], 1, 'omitnan');
max_max_density_over_sets = max(result_mat(:, :, 4), [], 1, 'omitnan');
mean_max_density_over_sets = mean(result_mat(:, :, 4), 1, 'omitnan');
%}

for j=1:exp_num
    fprintf("(DEPZ-NR)For exp: %d, with dir: %s , b: %.6f has time: %.1f seconds and max memory over cases: %.2e\n\n", ...
        j, mat2str(dirs{j}), bs{j}, total_time_per_exp(j), max_memory_per_exp(j));
    %{
    fprintf("Density info: min_init_d: %d, max_init_d: %d, mean_init_d: %.4f," + ...
        " min_max_d: %d, max_max_d: %d, mean_max_d: %.4f\n\n",...
        min_init_density_over_sets(j), max_init_density_over_sets(j), mean_init_density_over_sets(j),...
        min_max_density_over_sets(j), max_max_density_over_sets(j), mean_max_density_over_sets(j));
    %}
end


function [G, E, center, GI, alpha_size] = preprocess_data(file_path_E, file_path_G, file_path_c, file_path_GI, power_case)
    % Load and process matrix E.
    tmp = load(file_path_E);
    E = tmp.E';
    E = power_case * E;
    
    % Load and process matrix G.
    tmp = load(file_path_G);
    G = tmp.G';
    
    % Load the center.
    tmp = load(file_path_c);
    center = tmp.c;
    
    % Load and process matrix GI.
    tmp = load(file_path_GI);
    GI = tmp.GI';
    
    % Set alpha_size as the number of columns of E.
    alpha_size = size(E, 2);
end


function [a, b, adjusted_vector, adjusted_value, c, b_val, max_depth] = prepare_problem(exp_info, alpha_size, GI, center)
    % Extract experiment parameters.
    c = exp_info{1};      % Direction vector (assumed 1×p)
    b_val = exp_info{2};  % Scalar value
    max_depth = exp_info{3};

    % Define the interval bounds.
    a = -1*ones(1, alpha_size);
    b = ones(1, alpha_size);

    % Preprocess the squared terms into positive domain only TBD

    % Prepare the adjusted vector and value.
    % GI is assumed to be m×alpha_size and c is 1×p (with appropriate dimensions).
    % Here we compute: alpha = -sign(GI * c') and then
    % alpha_sum_vector = sum(alpha .* GI, 1)
    alpha = -sign(GI * c');  
    alpha_sum_vector = sum(alpha .* GI, 1);
    adjusted_vector = alpha_sum_vector + center(:)';  % Ensure center is a row vector.
    adjusted_value = dot(adjusted_vector, c);
end


function [beta_min, beta_max] = overapprox_depz(G, a, b, E)
    % a, b are row vectors 
    a_exp = a .^ E;  
    b_exp = b .^ E;

    % Common case for the max should choose the bigger one.
    beta_max = max(a_exp(:, 1), b_exp(:,1));

    % For the beta_min need to distinguish the case the domain crossing 0
    beta_min = min(a_exp(:, 1), b_exp(:,1));
    if a(1) * b(1) < 0
        assert (a(1) < 0 && b(1) > 0)
        % Get the cases where a_exp is bigger than zero. Note we should
        % also take care of the case where -1^0 == 1 is applied
        reset_indices = find(a_exp(:, 1) > 0);
        beta_min(reset_indices) = 0;

        % For those not existing indices, need to resect back to 1
        not_exit_indices = find(E(:, 1) == 0);
        beta_min(not_exit_indices) = 1;
    end

    for i = 2:size(a_exp,2)
        temp_max = max(a_exp(:, i), b_exp(:, i));
        temp_min = min(a_exp(:, i), b_exp(:, i));
        if a(i)*b(i) < 0
            temp_reset_indices = find(a_exp(:, i) > 0);
            temp_min(temp_reset_indices) = 0;

            temp_not_exit_indices = find(E(:, i) == 0);
            temp_min(temp_not_exit_indices) = 1;
        end

        candidates = [temp_min .* beta_min, temp_max .* beta_min, ...
                      temp_min .* beta_max, temp_max .* beta_max];

        beta_min = min(candidates, [], 2);
        beta_max = max(candidates, [], 2);
    end
end


function point = middle_point_polynomial_zonotope_with_dom(G, a, b, E, adjusted_vector, dir)
    % Compute the middle coefficients.
    middle_coefficients = (a + b) / 2;
    
    % Compute the transformed coefficients by taking the product over each row.
    left_transformed   = prod(a .^ E, 2);
    middle_transformed = prod(middle_coefficients .^ E, 2);
    right_transformed  = prod(b .^ E, 2);
    
    % Compute the zonotope points. Here we assume that G has size m×p, so that
    % left_transformed' * G produces a 1×p row vector.
    left_point   = left_transformed' * G + adjusted_vector;
    middle_point = middle_transformed' * G + adjusted_vector;
    right_point  = right_transformed' * G + adjusted_vector;
    
    % Compute inner products with the given direction vector.
    left_inner  = left_point * dir';
    middle_inner = middle_point * dir';
    right_inner = right_point * dir';
    
    % Choose the point with the smallest inner product.
    values = [left_inner, middle_inner, right_inner];
    [~, idx] = min(values);
    if idx == 1
        point = left_point;
    elseif idx == 2
        point = middle_point;
    else
        point = right_point;
    end
end


function result = depz_intersection(G, a, b, E, c, b_val, adjusted_vector, adjusted_value, mem_track, depth, max_depth, density_info)
    
    % Terminate if maximum recursion depth is exceeded.
    if depth > max_depth
        result = [];  % Return empty to indicate an undecidable case.
        return;
    end
    
    % Compute interval over-approximations.
    [beta_min, beta_max] = overapprox_depz(G, a, b, E);
    
    % Check the half-space intersection.
    [intersect, intersect_val, biggest_gen_intv_idx] = check_half_space_intersection(E, G, beta_min, beta_max, c, b_val, adjusted_value);
    
    if ~intersect
        result = false;
        return;
    end
    
    % Compute a candidate (middle) point.
    mid_pt = middle_point_polynomial_zonotope_with_dom(G, a, b, E, adjusted_vector, c);
    if dot(c, mid_pt) <= b_val
        result = true;
        return;
    end
    
    
    % A heuristic based method of determining the factor to split, based on
    % interval range. Note as long as multiple factor exit in the geneartor
    % term, after one is splitted, another would be chosen based on the
    % largest range. 
    non_zero_factor_indices = find(E(biggest_gen_intv_idx, :) > 0);
    t_low = zeros(size(E, 2));
    t_low(non_zero_factor_indices) = a(non_zero_factor_indices);

    t_up = zeros(size(E, 2));
    t_up(non_zero_factor_indices) = b(non_zero_factor_indices);
    
    [~, split_index] = max(t_up - t_low);
    

    % Perform the cyclic splitting.
    split_a_1 = a;
    split_b_1 = b;
    split_b_1(split_index) = (a(split_index) + b(split_index)) / 2;
    
    split_a_2 = a;
    split_b_2 = b;
    split_a_2(split_index) = (a(split_index) + b(split_index)) / 2;
    
    % Update the memory tracker and density (though density will remain constant, still tracking).
    new_mem_usage = (numel(split_a_1) + numel(split_b_1)) + (numel(split_a_2) + numel(split_b_2));
    mem_track.add(new_mem_usage);
    
    % Update the density info
    [temp_density, ~] = size(E);
    if temp_density > density_info.get_max()
        density_info.update_max(temp_density);
        fprintf("New max density is %d\n", density_info.get_max());
    end
    
    % Recursively check the two splits.
    result_1 = depz_intersection(G, split_a_1, split_b_1, E, c, b_val, adjusted_vector, adjusted_value, ...
                                 mem_track, depth + 1, max_depth, density_info);
    if (isequal(result_1, true) || isempty(result_1))
        result = result_1;
        return;
    end
    
    result_2 = depz_intersection(G, split_a_2, split_b_2, E, c, b_val, adjusted_vector, adjusted_value, ...
                                 mem_track, depth + 1, max_depth, density_info);
    
    if (isequal(result_2, true) || isempty(result_2))
        result = result_2;
        return;
    end
    
    assert (isequal(result_1, false) && isequal(result_2, false));
    result = false;
    return;
end


function [intersect, val, biggest_gen_intv_idx] = check_half_space_intersection(E, G, beta_min, beta_max, c, b, adjusted_value)
    % Compute the product G*c.
    cG = G * c(:);  % Ensure c is a column vector.
    
    % For each element of cG, if it is nonnegative use beta_min, else beta_max.
    temp = zeros(size(cG));
    
    % Get the positive generators, use beta_min to minimize the sum
    pos_idx = (cG >=0);

    % Get the negative generators, use beta_max to minimize the sum
    neg_idx = (cG<0);

    % Update the corresponding evaluated G matrices
    temp(pos_idx) = beta_min(pos_idx) .* cG(pos_idx);
    temp(neg_idx) = beta_max(neg_idx) .* cG(neg_idx);
    min_value = sum(temp);
    
    % Check the half-space condition.
    val = min_value + adjusted_value;
    intersect = (val <= b);

    % Additional function telling the biggest multiplied interval
    lower_intv = zeros(size(cG));
    upper_intv = zeros(size(cG));

    lower_intv(pos_idx) = beta_min(pos_idx) .* cG(pos_idx);
    upper_intv(pos_idx) = beta_max(pos_idx) .* cG(pos_idx);

    lower_intv(neg_idx) = beta_max(neg_idx) .* cG(neg_idx);
    upper_intv(neg_idx) = beta_min(neg_idx) .* cG(neg_idx);
   
    assert(all(lower_intv < upper_intv));
    [~, biggest_gen_intv_idx] = max(upper_intv - lower_intv); 
end

function print_result(check_result, i, j, time_usage, mem_track)
    time_str = sprintf(' Time usage is: %fs. with memory usage: %d', time_usage, mem_track.get());
    if isequal(check_result, true)
        fprintf('Case %d, Exp %d: Intersection found with the hyperplane.%s\n', i, j, time_str);
        error_mesg = sprintf('Case %d, Exp %d: Intersection found with the hyperplane.%s\n', i, j, time_str);
        assert(false, error_mesg);
    elseif isequal(check_result, false)
        fprintf('Case %d, Exp %d: No intersection found.%s\n', i, j, time_str);
    elseif isempty(check_result)
        fprintf('Case %d, Exp %d: Inconclusive result - reached maximum recursion depth.%s\n', i, j, time_str);
        error_mesg = sprintf('Case %d, Exp %d: Inconclusive result - reached maximum recursion depth.%s\n', i, j, time_str);
        assert(false, error_mesg);
    end
end