function mtaff_new_plot(params)
% ========================================================================
%{
Input: elements stored inside a params structure
    given_pZ: the input pz, in the case of var_dubins would be the veri pz
    pzfilename: the filename composing the desired plotted pz
    dirfilename: the filename storing the directions and their optimal pts
    origin_splits_num: number of splits for the default plotting technique
    itct_split_num: number of splits used for the intersection with new
                    method
    projdim: the plotting project dimension, since there are 2 systems
             though both of them working on project on the first dim
    option: "default" -- plot the PZ though both default + random points
                     Record: default computation time + new method time
            
            "countor" -- plot the countor of default + new method + 
                         exact extremes + intersection of supports 
                     Record: since demonstration purpose, no record of time
            
            "test" -- plot the default + new method, 
                     Record: no record of time
    rdpts_num: number of random points to generate for the "default"
    figfname: filename of the existing figure to compare with our method
              under the option for "compare"
    
%}   
% ========================================================================
% Load the parameters:
given_pZ = params.given_pZ;
pzfilename = params.pzfilename;
dirfilename = params.dirfilename;
origin_splits_num = params.origin_splits_num;
itct_split_num = params.itct_split_num;
projdim = params.projdim;
option = params.option;

% Follow parameters only for default version
if strcmp(option, "default") 
    rdpts_num = params.rdpts_num;
end

if strcmp(option, "compare")
    figfname = params.figfname;
end
% =========================================================================

% Shared Computation Part =================================================
% This part should not have any heavy computation
warning('off','all');
load(pzfilename);
load(dirfilename);
pZ_set = polyZonotope(pz.c, pz.G, pz.Grest, pz.expMat, pz.id);
assert(isequal(given_pZ, pZ_set), ['The calculate pZ does not match with' ...
    ' the offline stored pZ for extreme points calculation']);

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

% =========================================================================
% Below are the different visualization based on the options

if strcmp(option, "compare")
    % Get the existing figure to compare with
    openfig(figfname, 'visible');
    hold on;

    % Get the intersected polygon by our method
    diaryname = "multiaff/out/dubins/compare_exp_log.txt";
    compare_figname = "multiaff/out/dubins/compare_figure.fig";
    diary(diaryname);
    
    % Compute and plot our method for comparison
    starttime = datetime('now');
    disp(['Start time of computing + plot intersection usage pz with splitnum of: ' ...
        ,num2str(itct_split_num) ,' splits is: ', char(starttime)]);
    
    % Main block of computing the intersection
    
    [~, xVals, yVals] = new_pgon_pts(pZ_set, projdim, 'EdgeColor','none', ...
        'Splits',itct_split_num);
    default_pz_pgon = polyshape(xVals, yVals);
    extreme_pz_pgon = polyshape(complete_contour_pts(1, complete_ch_idx),...
        complete_contour_pts(2, complete_ch_idx));
    intersect_pz_pgon = intersect(default_pz_pgon, extreme_pz_pgon);
    plot(intersect_pz_pgon, 'FaceColor', 'none', 'FaceAlpha',1, ...
        'EdgeColor','b');
    %legend(hanour, "Our Method");
    saveas(gcf, compare_figname);
    
    % Output the ending time
    endtime = datetime('now');
    disp(['End time of computing + plot intersection usage pz with splitnum of: ' ...
        ,num2str(itct_split_num) ,' splits is: ',  char(endtime)]);
    
    % Close the log file
    diary off;
end

if strcmp(option, "test")
    figure; hold on;
    plot(pZ_set, projdim,'FaceColor','none', 'EdgeColor','r', ...
        'Splits',origin_splits_num);
    [~, xVals, yVals] = new_pgon_pts(pZ_set, projdim, 'EdgeColor','none', ...
        'Splits',itct_split_num);
    default_pz_pgon = polyshape(xVals, yVals);
    extreme_pz_pgon = polyshape(complete_contour_pts(1, complete_ch_idx),...
        complete_contour_pts(2, complete_ch_idx));
    intersect_pz_pgon = intersect(default_pz_pgon, extreme_pz_pgon);
    plot(intersect_pz_pgon, 'FaceColor', 'none', 'FaceAlpha',1, ...
        'EdgeColor','b');
    xlim([-15, 15]);
    ylim([4, 16]);
    xlabel(['x_{',num2str(1),'}']);
    ylabel(['x_{',num2str(2),'}']);
    set(gcf, 'Units', 'centimeters', 'OuterPosition', [0, 0, 14, 10]);
    box;
end    

if strcmp(option, "default")
    ptsnum = rdpts_num;
    splitnum = origin_splits_num;
    diaryname = sprintf('%s%d%s%d%s', "multiaff/out/dubins/default_plot_text_split_", ...
        splitnum, "_rdpts_",ptsnum, '.txt');
    figname = sprintf('%s%d%s%d%s', "multiaff/out/dubins/default_plot_split_", ...
        splitnum, "_rdpts_",ptsnum, '.fig');
    
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
    legend2 = sprintf('%s', 'PZ Random Points');
    legend([hanPZ, hanpts], legend1, legend2);
    set(gcf, 'Units', 'centimeters', 'OuterPosition', [0, 0, 14, 10]);
    box;
    saveas(gcf, figname);
    diary off;
end

if strcmp(option, "countor")
    figure; hold on;
    [~, xVals, yVals] = new_pgon_pts(pZ_set, [1 2],'Splits',origin_splits_num);
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
end

end

