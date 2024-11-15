% Initialize the folder path containing the .mat files
folderPath = 'laubLoomis';

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
split = 3;
start_idx = 1;
end_idx = 2000;
exp_num = length(dirs);
result_mat = zeros(end_idx, exp_num, 2); % first store time, second store memory

% Iterate over all cases
for i = start_idx:end_idx
    % Load G and E data from the respective files
    G_data = load(fullfile(folderPath, sprintf('G_interval_%d.mat', i)));
    E_data = load(fullfile(folderPath, sprintf('E_interval_%d.mat', i)));
    c_data = load(fullfile(folderPath, sprintf('c_interval_%d.mat', i)));
    GI_data = load(fullfile(folderPath, sprintf('GI_interval_%d.mat', i)));
    % Extract G and E from loaded data (assuming they are stored with known variable names)
    G = G_data.G;  % Adjust if the variable name inside the file is different
    E = E_data.E;  % Adjust if the variable name inside the file is different
    c = c_data.c;
    GI = GI_data.GI;
    % Create the PolyZonotope object pZ
    pZ = polyZonotope(c, G, GI, 2.*E);
    init_mem = (size(G,1)*size(G, 2)) + (size(E,1)*size(E, 2)) + (size(GI,1)*size(GI, 2));
    for j = 1:exp_num
        dir = dirs{j};
        b = bs{j};
        % Create the halfspace object hs1
        hs1 = halfspace(dir, b);
        
        % Check intersection
        start_time = tic;
        [a, mem] = isIntersecting_(pZ, hs1, 'approx', split);
        tComp = toc(start_time);
        fprintf("Set %d, Exp %d, has time %s, memory %s\n", i, j, num2str(tComp), num2str(init_mem + mem));
        result_mat(i, j, 1) = tComp;
        result_mat(i, j, 2) = init_mem + mem;
    end
    
end

% Save the result_mat
save('cora_laub_mem_time.mat', 'result_mat');