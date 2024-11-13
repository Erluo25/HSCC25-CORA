function better_pZ_plot(given_pZ, given_zono, end_time_points, pzfilename, ...
    dirfilename, origin_splits_num, itct_split_num, projdim)
warning('off','all');

% Load the offline computed extreme points for the given_pZ
load(pzfilename);
load(dirfilename);
pZ_set = polyZonotope(pz.c, pz.G, pz.Grest, pz.expMat, pz.id);
%assert(isequal(given_pZ, pZ_set), ['The calculate pZ does not match with' ...
%    ' the offline stored pZ for extreme points calculation']);
%assert(isequal(round(given_zono.Z, 4), round(zono_set.Z, 4)), ['The calculate zono does not match with' ...
%    ' the offline stored zono']);


% Trying to plot the PZ by default together with the random points=========
%{
ptsnum = 1000000;
splitnum = 50;
diaryname = sprintf('%s%d%s%d%s', "default_plot_text_split_", splitnum,...
    "_rdpts_",ptsnum, '.txt');
figname = sprintf('%s%d%s%d%s', "default_plot_split_", splitnum,...
    "_rdpts_",ptsnum, '.fig');

diary(diaryname);
figure; hold on;
% Get the random points and plot them
tic
rdpts = randPoint(given_pZ, ptsnum);
tComp = toc;
disp(['computation time for getting ',num2str(ptsnum),' random point is: ',num2str(tComp), ' seconds']);

% Plot the PZ 
starttime = datetime('now');
disp(['Start time of plotting with ',num2str(splitnum) ,' splits is: ', char(starttime)]);

% Main function of performing the plot
hanPZ = plot(given_pZ, [1 2],'Filled',true, 'FaceColor',[.8 .8 .8], ...
    'EdgeColor','none', 'Splits',splitnum);

endtime = datetime('now');
disp(['End time of plotting with ',num2str(splitnum) ,'splits is: ', char(endtime)]);

% Plot the rdpts to avoid overlapping
hanpts = plot(rdpts(1, :),rdpts(2, :),'.k','MarkerSize',10);

% Resize and modify the legend
xlim([-15, 15]);
ylim([4, 16]);
xlabel(['x_{',num2str(1),'}']);
ylabel(['x_{',num2str(2),'}']);
legend1 = sprintf('%s%d%s', 'PZ With ', splitnum, ' Splits');
legend2 = sprintf('%d%s', ptsnum, ' Random Points');
legend([hanPZ, hanpts], legend1, legend2);
set(gcf, 'Units', 'centimeters', 'OuterPosition', [0, 0, 14, 10]);
box;
saveas(gcf, figname);
diary off;
%}
%==========================================================================

% Reshape the loaded extreme points and directions.
pts = pts';
dirs = dirs';

% Store the corresponding line information
slopes_and_biases = zeros(2, length(dirs));
aux_pts = zeros(length(pts(:, 1)), length(pts(1, :)));

for i = 1:length(dirs(1, :))
    % Compute the basic info of each line
    pt = pts(:, i);
    dir = dirs(:, i);

    % Need to rotate the normal vector by 90 degrees 
    dir = [dir(2, 1); -dir(1, 1)];
    slope = dir(2, 1) / dir(1, 1);
    bias = pt(2, 1) - slope * pt(1, 1);
    slopes_and_biases(:, i) = [slope; bias];
    
    % Compute the aux point for visualization purpose
    aux_pt = [1; slope + bias];
    aux_pts(:, i) = aux_pt;
end

% Now call the convex hull to sort the pts. 
ch_idx = convhull(pts');
pts = pts(:, ch_idx);
slopes_and_biases = slopes_and_biases(:, ch_idx);
aux_pts = aux_pts(:, ch_idx);

enclosing_pts = [];
% Any close pair, compute the intersection of there support line.
for i = 1:(length(pts(1, :))-1)
    slope1 = slopes_and_biases(1, i);
    bias1 = slopes_and_biases(2, i);
    slope2 = slopes_and_biases(1, i+1);
    bias2 = slopes_and_biases(2, i+1);
    
    pt_x = (bias2 - bias1) / (slope1 - slope2);
    pt_y = slope1 * pt_x + bias1;
    enclosing_pt = [pt_x; pt_y];
    enclosing_pts = [enclosing_pts, enclosing_pt];
end

% obtain the all the countour points
complete_contour_pts = [pts, enclosing_pts];
complete_ch_idx = convhull(complete_contour_pts');

% Get the curve for the constant parameter case
%const_x_s = [end_time_points(1,:)];
%const_y_s = [end_time_points(2,:)];

%}
% Visualization -----------------------------------------------------------
%{
% Fig 1: Zonotope + default PZ plot + exact extremes of PZ
figure; hold on;
hanZono = plot(zono_set,[1 2],'Filled',true,'FaceColor',[.6 .6 .6],...
    'EdgeColor','none');
hanPZ = plot(given_pZ, [1 2],'Filled',true, 'FaceColor',[.8 .8 .8], 'Splits',20);
hanExtreme = plot(pts(1,:), pts(2,:), '.r', 'MarkerSize', 10);
% Plot the tangent lines
%for j = 1:length(pts(1, :))
%    line = horzcat(pts(:, j), aux_pts(:, j));
%    plot(line(1, :), line(2, :));
%end
% Align figure
xlim([-15, 15]);
ylim([4, 16]);
xlabel(['x_{',num2str(1),'}']);
ylabel(['x_{',num2str(2),'}']);
legend([hanZono, hanPZ, hanExtreme], 'Zonotope Method', 'Default PZ',...
    'PZ Exact Extremes');
set(gcf, 'Units', 'centimeters', 'OuterPosition', [0, 0, 14, 10]);
box;

% Fig 2: Exact Extremes + tangent lines + counter pts
figure; hold on;
%hanPZ = plot(given_pZ, [1 2],'Filled',true, 'FaceColor',[.8 .8 .8], 'Splits',20);
hanExtremePTS = plot(pts(1,:), pts(2,:), '.r', 'MarkerSize', 10);
% Plot the tangent lines
for j = 1:length(pts(1, :))
    line = horzcat(pts(:, j), aux_pts(:, j));
    plot(line(1, :), line(2, :));
end
hanEnclose = plot(enclosing_pts(1, :), enclosing_pts(2, :), '.g', 'MarkerSize',10);
xlim([-15, 15]);
ylim([4, 16]);
xlabel(['x_{',num2str(1),'}']);
ylabel(['x_{',num2str(2),'}']);
legend([hanExtremePTS, hanEnclose], 'PZ Exact Extremes', 'Tangent Intersections');
set(gcf, 'Units', 'centimeters', 'OuterPosition', [0, 0, 14, 10]);
box;
%}

% Fig 3: Default Polygon + Extreme Polygon + Counter Points
%{
figure; hold on;
[~, xVals, yVals] = pZ_plot_pgon_pts(pZ_set, [1 2],'Filled', true,'FaceColor',...
    [.8 .8 .8], 'EdgeColor','none','Splits',20);
default_pz_pgon = polyshape(xVals, yVals);
extreme_pz_pgon = polyshape(complete_contour_pts(1, complete_ch_idx),...
    complete_contour_pts(2, complete_ch_idx));
hanDefaultPG = plot(default_pz_pgon, 'FaceColor', 'none', 'FaceAlpha',1, ...
    'EdgeColor','k');
hanExtremePG = plot(extreme_pz_pgon, 'FaceColor', 'none', 'FaceAlpha',1, ...
    'EdgeColor','b');
hanExtPTS =  plot(pts(1, :), pts(2, :), '.g', 'MarkerSize', 10);
hanEnclosePTS = plot(enclosing_pts(1, :), enclosing_pts(2, :), '.r', 'MarkerSize', 10);
xlim([-15, 15]);
ylim([4, 16]);
xlabel(['x_{',num2str(1),'}']);
ylabel(['x_{',num2str(2),'}']);
legend([hanDefaultPG, hanExtremePG, hanExtPTS, hanEnclosePTS],'PZ Polygon (Default)', 'Extreme Countor Polygon', ...
    'PZ Exact Extremes', 'Intersection Of Supports');
set(gcf, 'Units', 'centimeters', 'OuterPosition', [0, 0, 14, 10]);
box;
%}

% Original Plot ==========================================================
figure; hold on;

% Plot the zonotope set at the end time point for comparison.
%hanZono = plot(zono_set,[1 2],'Filled',true,'FaceColor',[.6 .6 .6],...
%    'EdgeColor','none');

% Obtain the polygon of the standard pZ plotting method and use its
% intersection with the pZ plot computed by extreme points enclosure method
% to have a more accurate polynomial zonotope plot
plot(pZ_set, projdim,'Filled', true,'FaceColor',...
    'none', 'EdgeColor','r','Splits',origin_splits_num)
[~, xVals, yVals] = pZ_plot_pgon_pts(pZ_set, projdim,'Filled', true,'FaceColor',...
    [.8 .8 .8], 'EdgeColor','none','Splits',itct_split_num);
default_pz_pgon = polyshape(xVals, yVals);
extreme_pz_pgon = polyshape(complete_contour_pts(1, complete_ch_idx),...
    complete_contour_pts(2, complete_ch_idx));
intersect_pz_pgon = intersect(default_pz_pgon, extreme_pz_pgon);
hanPZ = plot(intersect_pz_pgon, 'FaceColor', 'none', 'FaceAlpha',1, ...
    'EdgeColor','b');

% Plot the end point set with constant parameters
%hd_exact_tp = plot(const_x_s, const_y_s, 'r','LineWidth',2);

% Visualizing the extreme points and their corresponding support functions
%plot(pts(1,:), pts(2,:), '.r', 'MarkerSize', 10);
%plot(enclosing_pts(1, :), enclosing_pts(2, :), '.g', 'MarkerSize',10);
%plot(complete_contour_pts(1, complete_ch_idx), complete_contour_pts(2, complete_ch_idx), 'b');
%for j = 1:length(pts(1, :))
%    line = horzcat(pts(:, j), aux_pts(:, j));
%    plot(line(1, :), line(2, :));
%end
% ==========================================================================
%xlim([-15, 15]);
%ylim([4, 16]);
%xlabel(['x_{',num2str(1),'}']);
%ylabel(['x_{',num2str(2),'}']);
%legend([hanZono, hanPZ, hd_exact_tp], 'Zonotope Method', 'Our Approach',...
%    'Constant Parameter');
%legend([hanPZ, hd_exact_tp], 'Zonotope Method', 'Our Approach',...
%    'Constant Parameter');

%set(gcf, 'Units', 'centimeters', 'OuterPosition', [0, 0, 14, 10]);
box;
%}
end

