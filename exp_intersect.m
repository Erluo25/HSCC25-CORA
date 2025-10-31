function [i, a, result_mat] = exp_intersect(case_num, dir, b, splits, pZlist, diary_name, diary_folder)
warning('off', 'all');
if strlength(diary_name) > 0
    filepath = fullfile(diary_folder, diary_name);
    diary(filepath);
end

fprintf("dir: %s,  b: %.6f, splits: %d\n\n", mat2str(dir), b, splits)

result_mat = nan(case_num, 4);
cum_time = 0;
went_over_case_num = 0;

min_init_density = inf;
max_init_density = -inf;
cum_init_density = 0;

min_max_density = inf;
max_max_density = -inf;
cum_max_density = 0;

% Iterate over all cases
for i = 1:case_num
    pZ = pZlist{i};
    % Create the halfspace object hs1
    hs1 = halfspace(dir, b);
    
    % Create pseduo density recording for the default CORA intersection
    % checking as we are not reporting that number.
    density_info.init_density = 0;
    density_info.max_density = 0;
    
    % Check intersection
    tStart = tic;
    [a, mem] = isIntersecting_(pZ, hs1, 'approx', splits);
    %[a, mem, density_info] = improved_benchmark(pZ, hs1, splits);
    time_usage = toc(tStart);
    
    % Save the results to the result_mat
    result_mat(i, 1) = time_usage;
    result_mat(i, 2) = mem;
    result_mat(i, 3) = density_info.init_density;
    result_mat(i, 4) = density_info.max_density;

    % Format the counter(s) and print them out
    cum_time = cum_time + time_usage;
    went_over_case_num = went_over_case_num + 1;
    
    % Go over the init_density cases
    min_init_density = min(min_init_density, density_info.init_density);
    max_init_density = max(max_init_density, density_info.init_density);
    cum_init_density = cum_init_density + density_info.init_density;
    mean_init_density = cum_init_density / went_over_case_num;

    % Go over the max_density cases
    min_max_density = min(min_max_density, density_info.max_density);
    max_max_density = max(max_max_density, density_info.max_density);
    cum_max_density = cum_max_density + density_info.max_density;
    mean_max_density = cum_max_density / went_over_case_num;

    % Display the information
    if strlength(diary_name) > 0
        fprintf("For case: %d, went over cases num: %d, has the cum_time: %.4f,\n" + ...
            "min_init_density: %d, max_init_density: %d, mean_init_density: %.4f\n" + ...
            "min_max_density: %d, max_max_density: %d, mean_max_density: %.4f\n\n",...
            i, went_over_case_num, cum_time, ...
            min_init_density, max_init_density, mean_init_density,...
            min_max_density, max_max_density, mean_max_density);
    end

    % Break the loop if intersection happens
    if a == 1
        break;
    end
end

if strlength(diary_name) > 0
    diary off;
end
