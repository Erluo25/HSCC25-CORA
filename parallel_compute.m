function parallel_compute(diary_output_fn, timeout, power_case, dirs, bs, ...
    split_list, folderPath, case_num, worker_diary)
% Set the warning and the diary output filename
warning('off', 'all');
diary(diary_output_fn);
output_args_num = 3;
diary_folder = '/Users/ertailuo/Research/CORA-2024/ParallelResults';

% Start a parallel pool if not already started
if isempty(gcp('nocreate'))
    parpool;
end

% Check the input dimensions are proper
assert((length(dirs) == length(bs)) && (length(bs) == length(split_list)), ...
    'invalid input lengths')

% Set up the number of jobs and the pZ list
numJobs = length(dirs);
pZlist = cell(case_num, 1);

% Preload the sets to speed up the parallel computation
for i = 1:case_num
    % load the data from the saved files of the set representations
    G_data = load(fullfile(folderPath, sprintf('G_interval_%d.mat', i)));
    E_data = load(fullfile(folderPath, sprintf('E_interval_%d.mat', i)));
    c_data = load(fullfile(folderPath, sprintf('c_interval_%d.mat', i)));
    GI_data = load(fullfile(folderPath, sprintf('GI_interval_%d.mat', i)));
    
    % Extract G and E from loaded data (assuming they are stored with known variable names)
    G = G_data.G;  
    E = E_data.E; 
    c = c_data.c;
    GI = GI_data.GI;
    
    % Create the PolyZonotope object pZ
    pZ = polyZonotope(c, G, GI, power_case.*E);
    pZlist{i} = pZ;
end

% Pre-allocate storage for FevalFuture objects and start times
futures = parallel.FevalFuture.empty(numJobs, 0);
startTime = cell(numJobs, 1);

% Allocate the result cell and result_mat(storing sparsity info for each case):
resultCell = repmat({"No result"}, 1, numJobs); 

% Start all jobs with parfeval and record their start times
for i = 1:numJobs
    % Load all the data necessary for the experiment
    dir = dirs{i};
    b = bs{i};
    split = split_list{i};
    diary_name = "";
    if worker_diary
        diary_name = sprintf("%s_exp_%d_dir_%s_b_%.6f_case.txt", folderPath, i, mat2str(dir), b);
    end

    futures(i) = parfeval(@exp_intersect, output_args_num, ...
                            case_num,dir,b,split, pZlist, diary_name, diary_folder);
    startTime{i} = tic;  % Record start time id for each job
end

% Track completed jobs
isJobFinished = false(numJobs, 1);  

% Monitor jobs until all are complete
while ~all(isJobFinished)           
    % Loop until all jobs are finished
    for i = 1:numJobs
        if isJobFinished(i)
            continue;  % Skip completed jobs
        end
        
        % Check the elapsed time for the current job
        elapsedTime = toc(startTime{i});
        
        % Check if the job is finished
        if strcmp(futures(i).State, 'finished')
            % Fetch the result if completed
            [idx, a, result_mat] = fetchOutputs(futures(i));
            
            % Compute the more accurate timing by getting rid of the
            % overheads
            more_accurate_timing = sum(result_mat(:, 1), 1, 'omitnan');

            resultCell{i} = sprintf("Intersection checking completed for exp: %d, with rough time: %s seconds, " + ...
                "and accurate timing: %.1f seconds, with i: %d, with check result a: %d, " + ...
                "dir: %s , b: %.6f, split num: %d \n", ...
                i, num2str(elapsedTime), more_accurate_timing, idx, a, mat2str(dirs{i}), bs{i}, split_list{i});
            isJobFinished(i) = true;  % Mark this job as finished
            
        elseif elapsedTime > timeout
            % Cancel the job if it exceeds the timeout
            cancel(futures(i));
            resultCell{i} = sprintf("Intersection checking cancelled for exp: %d, " + ...
                "with rough time: %s seconds, dir: %s , b: %.6f, split num: %d\n", ...
                i, num2str(elapsedTime), mat2str(dirs{i}), bs{i}, split_list{i});
            isJobFinished(i) = true;  % Mark this job as finished
        end
    end
    pause(0.5);  % Pause briefly to avoid excessive CPU usage
end

% Output the results
for j=1:numJobs
    fprintf(resultCell{j});
end
diary off;

end