%dataset = 'laubLoomis';
dataset = 'VanDelPol';
start_idx = 1;
end_idx = 1348;
% Try out unsquared depz
exp_info = {[1, 0], -2.026, 40};
exp_info = {[-1, 0], -2.143, 40};
exp_info = {[0, 1], -2.71, 40};
exp_info = {[0, -1], -2.79, 40};

% Experiment for depz ==================================
%exp_info = {[1, 0], -2.0165, 40};
%exp_info = {[-1, 0], -2.138, 40};
%exp_info = {[0, 1], -2.73, 40};
%exp_info = {[0, -1], -2.804, 40};
% Laub stuff========================
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
k = 6;
exp_info = {dirs{k}, bs{k}, 40};\
%}
% Laub stuff========================


cum_time = 0;
for i = start_idx:end_idx
    file_path_E  = fullfile(dataset, sprintf('E_interval_%d.mat', i));
    file_path_G  = fullfile(dataset, sprintf('G_interval_%d.mat', i));
    file_path_c  = fullfile(dataset, sprintf('c_interval_%d.mat', i));
    file_path_GI = fullfile(dataset, sprintf('GI_interval_%d.mat', i));
    [G, E, center, GI, alpha_size] = preprocess_data(file_path_E, file_path_G, file_path_c, file_path_GI);
    
    % Make up the direction and initialzied info
    [a, b, adjusted_vector, adjusted_value, c, b_val, max_depth] = prepare_problem(exp_info, alpha_size, GI, center);
    init_mem = numel(G) + numel(E) + numel(GI);
    mem_track.val = init_mem;
    
    % Time the intersection checking.
    tStart = tic;
    check_result = depz_intersection(G, a, b, E, c, b_val, adjusted_vector, adjusted_value, ...
                                     mem_track, 0, max_depth, 1, true, 1e-7);
    time_usage = toc(tStart);
    cum_time = cum_time + time_usage;
    print_result(check_result, i, time_usage);
    %disp(["Memory usage for depz is: %s", num2str(mem_track.val)]);
end
% Note print_result will panic if there's inconclusive result or intersect
fprintf("DEPZ time is: %d", cum_time);


% Experiment for CORA PZ =================================
%{
hs = halfspace([-1, 0], -2.14);
split_num = 6;
cum_time = 0;
for i = start_idx:end_idx
    file_path_E  = fullfile(dataset, sprintf('E_interval_%d.mat', i));
    file_path_G  = fullfile(dataset, sprintf('G_interval_%d.mat', i));
    file_path_c  = fullfile(dataset, sprintf('c_interval_%d.mat', i));
    file_path_GI = fullfile(dataset, sprintf('GI_interval_%d.mat', i));
    [G, E, center, GI, alpha_size] = preprocess_data(file_path_E, file_path_G, file_path_c, file_path_GI);
    pZ = polyZonotope(center,G', GI', E');
    
    % Check intersection
    start_time = tic;
    [cora_result, mem] = isIntersecting_(pZ, hs, 'approx', split_num);
    tComp = toc(start_time);
    cum_time = cum_time + tComp;
    fprintf("CORA: Set %d, Exp %d, has result %d time %s, memory %s\n", i, 1,...
        cora_result, num2str(tComp), num2str(init_mem + mem));
    if cora_result == 1
        fprintf("CORA: Set %d has intersection", i);
        break;
    end
end
fprintf("CORA total time: %d, split num: %d", cum_time, split_num);
% Plot the PZ =================================
%{
figure;
hold on;
plot(pZ, [1, 2], 'r', 'Splits',20);
%}
%}

function [G, E, center, GI, alpha_size] = preprocess_data(file_path_E, file_path_G, file_path_c, file_path_GI)
    % Load and process matrix E.
    tmp = load(file_path_E);
    E = tmp.E';
    %E = 2 .* E; 
    %{
    [~, factor_num] = size(E);
    selected_factor_num = min(factor_num, 10);
    E = 31 * E(:, 1:selected_factor_num);
    
    %E = 30 * E;
    %E = E(:, 1:6);
    %}
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
    %a = zeros(1, alpha_size);
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

%{
function [beta_min, beta_max] = overapprox_depz(G, a, b, E)
    % a and b are assumed to be 1×n row vectors.
    % E is an m×n matrix so that a.^E broadcasts to an m×n matrix.
    a_exp = a .^ E;  % Each row: a_i raised to the powers in E.
    b_exp = b .^ E;
    beta_min = prod(a_exp, 2);  % Product along each row.
    beta_max = prod(b_exp, 2);
end
%}

%{
function [beta_min, beta_max] = overapprox_depz_old(G, a, b, E)
    % a, b are row vectors 
    %disp('hi')
    a_exp = a .^ E;  
    b_exp = b .^ E;
    beta_min = min(a_exp(:,1), b_exp(:,1));
    beta_max = max(a_exp(:,1), b_exp(:,1));
    if a(1)*b(1) < 0
        beta_min = min(beta_min,0);
        beta_max = max(beta_max,0);
    end
    for i = 2:size(a_exp,2)
        candidates = [a_exp(:,i) .* beta_min, a_exp(:,i) .* beta_max, ...
                      b_exp(:,i) .* beta_min, b_exp(:,i) .* beta_max];
        beta_min = min(candidates, [], 2);
        beta_max = max(candidates, [], 2);
        if a(i)*b(i) < 0
            beta_min = min(beta_min,0);
            beta_max = max(beta_max,0);
        end
    end
end
%}

function [beta_min, beta_max] = overapprox_depz(G, a, b, E)
    % a, b are row vectors 
    %disp('hi')
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


function result = depz_intersection(G, a, b, E, c, b_val, adjusted_vector, adjusted_value, mem_track, depth, max_depth, split_index, repeat, tolerance)
    
    % Terminate if maximum recursion depth is exceeded.
    if depth > max_depth
        result = [];  % Return empty to indicate an undecidable case.
        return;
    end
    
    % Compute interval over-approximations.
    [beta_min, beta_max] = overapprox_depz(G, a, b, E);
    
    % Check the half-space intersection.
    [intersect, intersect_val, biggest_gen_intv_idx] = check_half_space_intersection(E, G, beta_min, beta_max, c, b_val, adjusted_value);
    %[intersect, intersect_val, ~] = check_half_space_intersection(G, beta_min, beta_max, c, b_val, adjusted_value);
    
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

    %[~, split_index] = max(b(non_zero_factor_indices) - a(non_zero_factor_indices));
    [~, split_index] = max(t_up - t_low);
    
    %{
    % Maybe we can try later splitting not based on b - a, but base on the
    % powered version
    temp_up = b.^E(biggest_gen_intv_idx, :);
    temp_low = a.^E(biggest_gen_intv_idx, :);
    up = max(temp_up, temp_low);
    low = min(temp_low, temp_up);

    cross_indices = find(a.*b < 0);
    low(cross_indices) = 0;
    
    zero_indices = find(E(biggest_gen_intv_idx, :)==0);
    low(zero_indices) = 1;

    [~, split_index] = max(up(non_zero_factor_indices) - low(non_zero_factor_indices));
    %}
    %{
    % Determine the index to split.
    % (Python code uses: split_index = (last_index+1) % length(a) with 0-based indexing.
    % Here we simulate that by computing a “python index” then converting to MATLAB’s 1-based index.)
    n = length(a);
    split_index = mod(split_index, n+1);
    if isequal(split_index, 0)
        split_index = 1;
    end
    %}

    % Perform the cyclic splitting.
    split_a_1 = a;
    split_b_1 = b;
    split_b_1(split_index) = (a(split_index) + b(split_index)) / 2;
    
    split_a_2 = a;
    split_b_2 = b;
    split_a_2(split_index) = (a(split_index) + b(split_index)) / 2;
    
    % Update the memory tracker.
    mem_track.val = mem_track.val + (numel(split_a_1) + numel(split_b_1)) + ...
                                 (numel(split_a_2) + numel(split_b_2));
    %{
    if ~repeat
        split_index = split_index + 1;
        repeat = true;
    else
        split_index = split_index;
        repeat = false;
    end
    %}
    %split_index = split_index + 1;

    % Recursively check the two splits.
    % Note: We pass last_index as a “python index” (i.e. MATLAB index-1).
    result_1 = depz_intersection(G, split_a_1, split_b_1, E, c, b_val, adjusted_vector, adjusted_value, ...
                                 mem_track, depth + 1, max_depth, split_index, repeat, tolerance);
    if isequal(result_1, true)
        result = true;
        return;
    end

    if isempty(result_1)
        result = [];
        return;
    end
    
    result_2 = depz_intersection(G, split_a_2, split_b_2, E, c, b_val, adjusted_vector, adjusted_value, ...
                                 mem_track, depth + 1, max_depth, split_index, repeat, tolerance);
    
    if isequal(result_2, true)
        result = true;
        return;
    end

    if isempty(result_2)
        result = [];
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
   
    %if(~all(lower_intv < upper_intv))
    %    assert(all(lower_intv < upper_intv));
    %end
    assert(all(lower_intv < upper_intv));
    [~, biggest_gen_intv_idx] = max(upper_intv - lower_intv); 
end

function print_result(check_result, i, time_usage)
    time_str = sprintf(' Time usage is: %fs.', time_usage);
    if isequal(check_result, true)
        fprintf('Case %d: Intersection found with the hyperplane.%s\n', i, time_str);
        error_mesg = sprintf('Case %d: Intersection found with the hyperplane.%s\n', i, time_str);
        assert(false, error_mesg);
    elseif isequal(check_result, false)
        fprintf('Case %d: No intersection found.%s\n', i, time_str);
    elseif isempty(check_result)
        fprintf('Case %d: Inconclusive result - reached maximum recursion depth.%s\n', i, time_str);
        error_mesg = sprintf('Case %d: Inconclusive result - reached maximum recursion depth.%s\n', i, time_str);
        assert(false, error_mesg);
    end
end