% Produce the non-monotonic decreasing figure
pZ = polyZonotope([0;0], [-2 -1 2; 0 -2 2], [0;0], [2 0 2; 0 2 6]);

% Define the split values to iterate over (one can try [6, 7] as well)
split_values = [0, 1, 2, 3];

% Define the ground truth split value
ground_truth_split = 13;

% Define a colormap for different splits
colors = lines(length(split_values));
ground_truth_fill_color = 'b'; 

% Create a new figure
figure;

% Hold on to plot all splits in the same figure
hold on;

% Loop over each split value to plot
for idx = 1:length(split_values)
    i = split_values(idx);
    % Assign a unique color for each split from the colormap
    current_color = colors(idx, :);
    % Plot the current split zonotope
    % Assuming pZ.plot accepts dimensions, color, and 'Splits' as a name-value pair
    pZ.plot([1, 2], 'Color', current_color, 'Splits', i);
end

% Plot the ground truth zonotope and fill it with blue
pZ.plot([1, 2], ground_truth_fill_color, 'Splits', ground_truth_split, 'FaceAlpha', 0.3);

% Add title to the figure
title('Polynomial Zonotope - All Splits with Unique Colors and Ground Truth');

% Label the axes
xlabel('Dimension 1');
ylabel('Dimension 2');

% Add legend to differentiate the plots
legend_strings = arrayfun(@(x) ['Splits = ' num2str(x)], split_values, 'UniformOutput', false);
legend_strings{end+1} = ['Ground Truth'];
legend(legend_strings, 'Location', 'best');

% Enhance grid for better visualization
grid on;

% Release the hold
hold off;