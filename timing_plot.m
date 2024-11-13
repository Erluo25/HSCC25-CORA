% Initialize the folder path containing the .mat files
folderPath = 'laubLoomis';
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
    pZ = polyZonotope(c, G, GI, 2.*E);
    plot(pZ, [1, 2], 'Splits', 8);
    % Define c1 and c2 satisfying the conditions
    % c1 = [1, 0];
    % c2 = [0, 1];
    
    % Define b (assuming a sample value for b; this should be provided based on your use case)
    %b = -2.138; % Modify as necessary for specific cases
    
    % Create the halfspace object hs1
    %hs1 = halfspace([-1, 0], b);
    
    % Check intersection
    %a = isIntersecting_(pZ, hs1, 'approx');
    %if a== 1
    %    break;
    %end
end

% Measure total elapsed time
total_time = toc(total_tic);

% Display the total elapsed time
fprintf('Total elapsed time for all cases: %.4f seconds\n', total_time);
