% run_experiment_laubLoomis.m
% This script sets up and runs the experiments for the 'laubLoomis' dataset.

dataset = 'laubLoomis';

% Define experiments as a cell array of cells.
% Each experiment is given as {direction, b_value, max_depth}.
depz_exps = { { [0, 0, 0, 0, 0, 0, 1], 0.137, 40}, ...
              { [0, 0, 0, 0, 1, 0, 0], 0.0685, 40}, ...
              { [1, 0, 0, 0, 0, 0, 0], 0.457, 40}, ...
              { [0, -1, 0, 0, 0, 0, 0], -1.354, 40}, ...
              { [0, 0, -1, 0, 0, 0, 0], -1.6505, 40}, ...
              { [0, 0, 0, 0, 0, 1, 0], 0.045, 40} };

start_idx = 1;
end_idx = 2000;  

% Run the depz experiment.
exp_depz_intersect(start_idx, end_idx, depz_exps, dataset, true, false);