% run_experiment.m
% This script sets up the dataset and experiment parameters and then calls the
% exp_depz_intersect function.

dataset = 'VanDelPol';

% Create cell arrays for experiments. Each experiment is defined as a cell array:
% {direction, b_val, max_depth}
depz_exps = { { [1, 0], -2.0165, 40}, {[-1, 0], -2.138, 40}, { [0, 1], -2.73, 40}, { [0, -1], -2.804, 40} };

start_idx = 1;
end_idx = 1348;

% Call the experiment function for depz.
exp_depz_intersect(start_idx, end_idx, depz_exps, dataset, true, false);
