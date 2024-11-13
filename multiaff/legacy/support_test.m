function support_test()

%Try out the 5d system taken from [1] with large size of sets
% References:
%    [1] Althoff, M.; Le Guernic, C. & Krogh, B. H. "Reachable Set
%        Computation for Uncertain Time-Varying Linear Systems", HSCC 2011

% System and parameters settings ------------------------------------------

dim_x = 5;
time_step = 0.05;
round = 3;

% Initialize the A matrix and dynamics
Acenter = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
Arad{1} = [0.1 0.1 0 0 0; 0.1 0.1 0 0 0; 0 0 0.1 0.1 0; 0 0 0.1 0.1 0; 0 0 0 0 0.1];
matZ_A = matZonotope(Acenter,Arad, 2);
dyn_sys = linParamSys(matZ_A, 1, 'varParam'); 

% Some share conributes like initial state and uncertain inputs
init_c = ones(dim_x,1);
init_G = 0.1*eye(dim_x);
input_c = zeros(dim_x,1);
input_G = 0.1*eye(dim_x);

% Set up the zonotope parameters
zono_params.tFinal = round * time_step;
zono_params.R0 = zonotope([init_c, init_G]);   
zono_params.U = zonotope([input_c,input_G]); 

zono_options.timeStep = 0.05;
zono_options.taylorTerms = 4;
zono_options.intermediateTerms = 2;
zono_options.zonotopeOrder = 20;

% Set up the PZ parameters
pZ_options.time_step = time_step;
pZ_options.round = round; 
pZ_options.break_order_mz = 2;
pZ_options.break_order_intv = 4;
pZ_options.max_order = 4;
pZ_options.input_break_order = 3;
pZ_options.pZ_order = 500;
pZ_options.restructure_ratio = 0.01;
pZ_options.restructure_order = 500;
pZ_options.reduce_idx = 7;
pZ_options.reduce_order = 50;
pZ_options.compact_flag = 1;

init_pZ = polyZonotope(init_c, init_G, [], eye(dim_x), [3, 4, 5, 6, 7]);
input_pZ = polyZonotope(input_c, input_G, [], eye(dim_x),...
    [8, 9, 10, 11, 12]);
setmaxid(1);
pZ_params.init = init_pZ;
pZ_params.input = input_pZ;

% Reachability Analysis ---------------------------------------------------

% Run the pz method
tic
rSet = reachParam(dyn_sys, pZ_params, pZ_options);
tComp = toc;
disp(['computation time of polynomial zonotope method: ',num2str(tComp)]);

dist_rSet = cell(length(rSet),1);
for i = 1:length(rSet)
    dist_rSet{i} = release_time(rSet{i});
end
save("HSCC_Examples_Experiments/testMtAff/support_test_sets.mat", "dist_rSet");

%{
pz = struct(dist_rSet{3});
save("HSCC_Examples_Experiments/testMtAff/test2.mat", 'pz');

%end_pz = release_time(rSet{end});

starttime = datetime('now');
disp(['Start time of plotting is: ', char(starttime)]);
better_pZ_plot(polyZonotope([0;0]), zonotope([]), [],'test2_red.mat', 'verts_test2_red_matlab.mat', 20, 5, [1 2]);

endtime = datetime('now');
disp(['End time of plotting is: ', char(endtime)]);
%}


end