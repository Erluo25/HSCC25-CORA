% Direct the output to the output_file
diary('laubLoomis_exp_output.txt');

% Define the timeout in seconds
timeout = 3; % change this to your desired timeout duration

% Start a parallel pool if not already started
if isempty(gcp('nocreate'))
    parpool;
end

output_args_num = 2;

folderPath = 'laubLoomis';
case_num = 2000;
dirs = {
    [0, 0, 0, 0, 0, 0, 1]; [0, 0, 0, 0, 0, 0, 1]; [0, 0, 0, 0, 0, 0, 1];
    [0, 0, 0, 0, 1, 0, 0]; [0, 0, 0, 0, 1, 0, 0];
    [1, 0, 0, 0, 0, 0, 0]; [1, 0, 0, 0, 0, 0, 0];
    [0, -1, 0, 0, 0, 0, 0];
    [0, 0, -1, 0, 0, 0, 0]; [0, 0, -1, 0, 0, 0, 0];
    [0, 0, 0, 0, 0, 1, 0]; [0, 0, 0, 0, 0, 1, 0]; [0, 0, 0, 0, 0, 1, 0];
};
bs = {
    0.137; 0.137; 0.137;
    0.0685; 0.0685;
    0.457; 0.457;
    -1.354;
    -1.6505; -1.6505;
    0.045; 0.045; 0.045;
    };
split_list = {
    8; 9; 10;
    6; 7;
    5; 6;
    25;
    11; 12;
    6; 7; 8;
};
assert((length(dirs) == length(bs)) && (length(bs) == length(split_list)), ...
    'invalid input lengths')

for i = 1:length(dirs)
    % Load all the data necessary for the experiment
    dir = dirs{i};
    b = bs{i};
    split = split_list{i};

    % Start the asynchronous execution of the function
    f = parfeval(@exp_intersect, output_args_num, ...
        folderPath,case_num,dir,b,split);
  
    % Start timing
    tic;
    while true
        % Check if the function is done
        if strcmp(f.State, 'finished')
            % Fetch the result if completed
            [idx, a] = fetchOutputs(f);
            tComp = toc;
            disp(['Intersection checking completed for exp:', num2str(i), ...
                ' with time: ', num2str(tComp), ' seconds with i: ', ...
                num2str(idx), ' has intersect: ', num2str(a), ' dir: ', num2str(dir), ...
                ' b:', num2str(b), ' split num: ', num2str(split)]);
            break;
        elseif toc > timeout
            % Cancel if it exceeds the timeout
            cancel(f);
            disp(['Intersection checking cancelled for exp:', num2str(i), 'with dir: ', num2str(dir), ...
                ' b:', num2str(b), ' split num: ', num2str(split), ...
                ' times out with tolarance: ', num2str(timeout), ' seconds']);
            break;
        else
            % Pause for a short duration to avoid busy-waiting
            pause(1);
        end
    end
end

diary off;  