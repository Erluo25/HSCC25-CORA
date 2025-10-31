% Direct the output to the output_file
diary_output_fn = 'van_exp_output.txt';
worker_diary = false;%true;

% Set the power case
power_case = 2;
%power_case = 31;

% Define the timeout in seconds
timeout = 900; % change this to your desired timeout duration
fprintf("time out is: %d seconds\n", timeout);

% Set up the experiments squared case
%
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
split_list = {
   10;
   6;
   3;
   2;
};
%}

% Set up the experiments of 31 degree cases
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

split_list = {
   40;
   40;
   40;
   40;
};
%}

% Set up the basic parameters
folderPath = 'VanDelPol';
case_num = 1348;

% Call the parallel function wrapper
parallel_compute(diary_output_fn, timeout, power_case, dirs, bs, split_list,...
    folderPath, case_num, worker_diary);