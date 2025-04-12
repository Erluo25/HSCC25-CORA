% Initialize the folder path containing the .mat files
folderPath = 'VanDelPol';

figure; hold on;
% Start timing the entire process
total_tic = tic;

% Iterate over all cases
for i = 1:1348
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
    %pZ = polyZonotope(c, G, GI, 2.*E);
    [factor_num, ~] = size(E);
    selected_factor_num = min(factor_num, 10);
    E = 31 * E(1:selected_factor_num, :);
    pZ = polyZonotope(c, G, GI, E);
    plot(pZ, [1, 2],'b', 'Splits', 7);
end

% Measure total elapsed time
total_time = toc(total_tic);

% Display the total elapsed time
fprintf('Total elapsed time for all cases: %.4f seconds\n', total_time);

% Plot all the halfspaces and save the zoomed figures.
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
exp_num = length(dirs);

for i=1:exp_num
    hs = halfspace(dirs{i}, bs{i});
    %plot(hs, [1, 2], 'r');
end


