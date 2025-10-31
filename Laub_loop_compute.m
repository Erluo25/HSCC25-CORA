% Initialize the folder path containing the .mat files
folderPath = 'laubLoomis';
%{
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
split = 40;
power_case = 2;
%}

% 31 cases =======================
power_case = 31;
dirs = {
    [0, 0, -1, 0, 0, 0, 0];
    [0, 0, 0, 0, 0, 1, 0];
};
bs = {
    -1.6508;
    -0.0501;
};
split = 40;

%{
dirs = {
    [0, 0, 0, 0, 0, 0, 1];
};
bs = {
    0.105; 
};
split = 40;

dirs = { 
    [0, 0, 0, 0, 1, 0, 0];
};
bs = {
    0.06;
};
split = 40;

dirs = {
    [1, 0, 0, 0, 0, 0, 0];
};
bs = {
    0.44;
};
split = 40;


dirs = {
    [0, -1, 0, 0, 0, 0, 0];
};
bs = {
    -1.373;
};
split = 40;

%}
%{
dirs = {
    [0, 0, -1, 0, 0, 0, 0];
};
bs = {
    -1.6508;
};
split = 40;

dirs = {
    [0, 0, 0, 0, 0, 1, 0];
};
bs = {
    -0.0501;
};
split = 40;
%}

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
    pZ = polyZonotope(c, G, GI, power_case*E);
   
    for j = 1:exp_num
        dir = dirs{j};
        b = bs{j};
        % Create the halfspace object hs1
        hs1 = halfspace(dir, b);
        
        % Check intersection
        start_time = tic;
        %[a, mem] = isIntersecting_(pZ, hs1, 'approx', split);
        [a, mem] = improved_benchmark(pZ, hs1, split);
        tComp = toc(start_time);
        assert(isequal(a, 0));
        fprintf("Set %d, Exp %d, intersects %d, has time %s, memory %s\n", i, j, a, num2str(tComp), num2str(mem));
        result_mat(i, j, 1) = tComp;
        result_mat(i, j, 2) = mem;
    end
end

% Print out the results
total_time_per_exp = sum(result_mat(:, :, 1), 1);
max_memory_per_exp = max(result_mat(:, :, 2), [], 1);
for j=1:exp_num
    fprintf("For exp: %d, with dir: %s , b: %s has time: %.1f seconds and max memory over cases: %.2e\n", ...
        j, mat2str(dirs{j}), num2str(bs{j}), total_time_per_exp(j), max_memory_per_exp(j));
end

% Save the result_mat
%save('cora_laub_mem_time.mat', 'result_mat');