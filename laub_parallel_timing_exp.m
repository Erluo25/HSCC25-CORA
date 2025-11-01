% Direct the output to the output_file
diary_output_fn ='laubLoomis_exp_output.txt';
worker_diary = false;%true;

% Set the power case
power_case = 2;
%power_case = 31;

% Define the timeout in seconds
timeout = 900; % change this to your desired timeout duration
fprintf("time out is: %d seconds\n", timeout);

% Set up the experiments squared version
%{
dirs = {
    [0, 0, 0, 0, 0, 0, 1]; [0, 0, 0, 0, 0, 0, 1];
    [0, 0, 0, 0, 1, 0, 0]; [0, 0, 0, 0, 1, 0, 0];
    [1, 0, 0, 0, 0, 0, 0]; [1, 0, 0, 0, 0, 0, 0];
    [0, -1, 0, 0, 0, 0, 0]; [0, -1, 0, 0, 0, 0, 0];
    [0, 0, -1, 0, 0, 0, 0]; [0, 0, -1, 0, 0, 0, 0];
    [0, 0, 0, 0, 0, 1, 0]; [0, 0, 0, 0, 0, 1, 0];
};
bs = {
    0.137; 0.137;
    0.0685; 0.0685;
    0.457; 0.457;
    -1.354; -1.354;
    -1.6505; -1.6505;
    0.045; 0.045;
    };
split_list = {
    9; 10;
    5; 6;
    4; 5;
    24; 25;
    10; 11;
    6; 7;
};
%}

% Also the squared version but only those cases where depth matching paper
%
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
split_list = {
    10;
    6;
    5;
    25;
    11;
    7;
};
%}

% Set up the experiments for cases with power 31
%{
dirs = {
    %[0, 0, 0, 0, 0, 0, 1]; 
    %[0, 0, 0, 0, 1, 0, 0];
    %[1, 0, 0, 0, 0, 0, 0];
    [0, -1, 0, 0, 0, 0, 0];
    %[0, 0, -1, 0, 0, 0, 0];
    %[0, 0, 0, 0, 0, 1, 0];
};
bs = {
    %0.105;
    %0.06;
    %0.44;
    -1.373;
    %-1.6508;
    %-0.0501;
};

split_list = {
    40;
    %40;
    %40;
    %40;
    %40;
    %40;
};
%}

% Set up the basic parameters
folderPath = 'laubLoomis';
case_num = 2000;

% Call the parallel function wrapper
parallel_compute(diary_output_fn, timeout, power_case, dirs, bs, split_list,...
    folderPath, case_num, worker_diary);