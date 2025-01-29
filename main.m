function main()

% Plot fig 3
figure3

saveas(gcf, 'Figure3.png')
close(gcf)

% Plot fig 4
timing_plot

% Save the original figure
saveas(gcf, 'Figure4.png');

% Zoom into the left
xlim([-2.05,-1.94]);
ylim([-0.2,0.2]);
saveas(gcf, 'Figure4-leftzoom.png');

% Zoom into the right
xlim([2.132, 2.15]);
ylim([-0.08, 0.1]);
saveas(gcf, 'Figure4-rightzoom.png');

close(gcf);

% Compute the table data for Van Del pol
van_parallel_timing_exp

% Compute the table data for Laubloomis
laub_parallel_timing_exp

% Compute the figure data for Van Del pol
Van_loop_compute

% Compute the figure data for Laubloomis
Laub_loop_compute


end