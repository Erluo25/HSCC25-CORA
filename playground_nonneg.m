dataset = 'VanDelPol';
i = 646;
file_path_E  = fullfile(dataset, sprintf('E_interval_%d.mat', i));
file_path_G  = fullfile(dataset, sprintf('G_interval_%d.mat', i));
file_path_c  = fullfile(dataset, sprintf('c_interval_%d.mat', i));
file_path_GI = fullfile(dataset, sprintf('GI_interval_%d.mat', i));
[G, E, center, GI, alpha_size] = preprocess_data(file_path_E, file_path_G, file_path_c, file_path_GI);

% Experiment for depz ==================================
% Make up the direction and initialzied info
exp_info = {[1, 0], -1.01, 40};
[a, b, adjusted_vector, adjusted_value, c, b_val, max_depth] = prepare_problem(exp_info, alpha_size, GI, center);
init_mem = numel(G) + numel(E) + numel(GI);
mem_track.val = init_mem;

% Time the intersection checking.
tStart = tic;
check_result = depz_intersection(G, a, b, E, c, b_val, adjusted_vector, adjusted_value, ...
                                 mem_track, 0, max_depth, 1, true, 1e-7);
time_usage = toc(tStart);

print_result(check_result, i, time_usage);
disp(["Memory usage for depz is: %s", num2str(mem_track.val)]);


% Experiment for CORA PZ =================================
pZ = polyZonotope(center,G', GI', E');
hs = halfspace([1, 0], -1.01);
% Check intersection
start_time = tic;
[cora_result, mem] = isIntersecting_(pZ, hs, 'approx', 10);
tComp = toc(start_time);
fprintf("CORA: Set %d, Exp %d, has result %d time %s, memory %s\n", 1, 1,...
    cora_result, num2str(tComp), num2str(init_mem + mem));

% Plot the PZ =================================
%{
figure;
hold on;
plot(pZ, [1, 2], 'r', 'Splits',20);
%}


function [G, E, center, GI, alpha_size] = preprocess_data(file_path_E, file_path_G, file_path_c, file_path_GI)
    % Load and process matrix E.
    tmp = load(file_path_E);
    E = tmp.E';
    E = 30 * E;
    E = E(:, 1:6);

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


function [beta_min, beta_max] = overapprox_depz(G, a, b, E)
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
    if ~check_half_space_intersection(G, beta_min, beta_max, c, b_val, adjusted_value)
        result = false;
        return;
    end
    
    % Compute a candidate (middle) point.
    mid_pt = middle_point_polynomial_zonotope_with_dom(G, a, b, E, adjusted_vector, c);
    if dot(c, mid_pt) <= b_val
        result = true;
        return;
    end
    
    % Determine the index to split.
    % (Python code uses: split_index = (last_index+1) % length(a) with 0-based indexing.
    % Here we simulate that by computing a “python index” then converting to MATLAB’s 1-based index.)
    n = length(a);
    split_index = mod(split_index, n+1);
    if isequal(split_index, 0)
        split_index = 1;
    end
    %python_index = mod(last_index + 1, n);  % Value in 0..(n-1)
    %split_index = python_index + 1;         % Convert to MATLAB index (1..n)
    
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
    split_index = split_index + 1;

    % Recursively check the two splits.
    % Note: We pass last_index as a “python index” (i.e. MATLAB index-1).
    result_1 = depz_intersection(G, split_a_1, split_b_1, E, c, b_val, adjusted_vector, adjusted_value, ...
                                 mem_track, depth + 1, max_depth, split_index, repeat, tolerance);
    if isequal(result_1, true)
        result = true;
        return;
    end
    result_2 = depz_intersection(G, split_a_2, split_b_2, E, c, b_val, adjusted_vector, adjusted_value, ...
                                 mem_track, depth + 1, max_depth, split_index, repeat, tolerance);
    
    if isequal(result_2, true)
        result = true;
        return;
    end

    % If either recursive call is undecidable (empty), return [].
    if isempty(result_1) || isempty(result_2)
        result = [];
    else
        result = false;
    end
end


function result = check_half_space_intersection(G, beta_min, beta_max, c, b, adjusted_value)
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
    result = (min_value + adjusted_value) <= b;
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