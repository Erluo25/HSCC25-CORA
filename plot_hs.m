hold on;
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
    plot(hs, [1, 2], 'r');
end


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
