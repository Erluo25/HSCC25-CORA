% Direct the output to the output_file
diary('laubLoomis_exp_output.txt');

% Define the timeout in seconds
timeout = 900; % change this to your desired timeout duration
fprintf("time out is: %d seconds\n", timeout);
% Start a parallel pool if not already started
if isempty(gcp('nocreate'))
    parpool;
end

% Set up the experiments
dirs = {
    [0, 0, 0, 0, 0, 0, 1]; [0, 0, 0, 0, 0, 0, 1];
    [0, 0, 0, 0, 1, 0, 0]; [0, 0, 0, 0, 1, 0, 0];
    [1, 0, 0, 0, 0, 0, 0]; [1, 0, 0, 0, 0, 0, 0];
    [0, -1, 0, 0, 0, 0, 0];
    [0, 0, -1, 0, 0, 0, 0];
    [0, 0, 0, 0, 0, 1, 0]; [0, 0, 0, 0, 0, 1, 0];
};
bs = {
    0.137; 0.137;
    0.0685; 0.0685;
    0.457; 0.457;
    -1.354;
    -1.6505;
    0.045; 0.045;
    };
split_list = {
    9; 10;
    6; 7;
    5; 6;
    25;
    11;
    7; 8;
};
assert((length(dirs) == length(bs)) && (length(bs) == length(split_list)), ...
    'invalid input lengths')

% Set up the basic parameters
output_args_num = 2;
folderPath = 'laubLoomis';
case_num = 2000;
numJobs = length(dirs);

% Pre-allocate storage for FevalFuture objects and start times
futures = parallel.FevalFuture.empty(numJobs, 0);
startTime = cell(numJobs, 1);

% Allocate the result cell:
resultCell = repmat({"No result"}, 1, numJobs); 

% Start all jobs with parfeval and record their start times
for i = 1:numJobs
    % Load all the data necessary for the experiment
    dir = dirs{i};
    b = bs{i};
    split = split_list{i};
    futures(i) = parfeval(@exp_intersect, output_args_num, ...
                            folderPath,case_num,dir,b,split);
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
            [idx, a] = fetchOutputs(futures(i));
            resultCell{i} = sprintf("Intersection checking completed for exp: %d, with time: %s seconds, " + ...
                "with i: %d, with check result a: %d, dir: %s , b: %s, split num: %d \n", ...
                i, num2str(elapsedTime), idx, a, mat2str(dirs{i}), num2str(bs{i}), split_list{i});
            isJobFinished(i) = true;  % Mark this job as finished
            
        elseif elapsedTime > timeout
            % Cancel the job if it exceeds the timeout
            cancel(futures(i));
            resultCell{i} = sprintf("Intersection checking cancelled for exp: %d, " + ...
                "with time: %s seconds, dir: %s , b: %s, split num: %d\n", ...
                i, num2str(elapsedTime), mat2str(dirs{i}), num2str(bs{i}), split_list{i});
            isJobFinished(i) = true;  % Mark this job as finished
        end
    end
    pause(1);  % Pause briefly to avoid excessive CPU usage
end

% Output the results
for j=1:numJobs
    fprintf(resultCell{j});
end
diary off;  