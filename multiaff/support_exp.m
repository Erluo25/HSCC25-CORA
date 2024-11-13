function support_exp(split_num, startidx)

% Load the common sets ====================================================
load('/Users/ertailuo/Research/CORA-2024/multiaff/data/supports/struct_supports.mat');
load('/Users/ertailuo/Research/CORA-2024/multiaff/data/supports/support_redset.mat');
% pZ_set comes from the struct of pz
pZ_set = polyZonotope(pz.c, pz.G, pz.Grest, pz.expMat, pz.id);
% given_pZ compes from the reduced pz
given_pZ = polyZonotope(pZred.c, pZred.G, pZred.Grest, pZred.expMat, pZred.id);
assert (isequal(pZ_set, given_pZ), 'Not equal sets');
% Project the pz on to [1 2] dimension
pZ_set = project(pZ_set, [1, 2]);
given_pZ = project(given_pZ, [1, 2]);

% Load the directions =====================================================
load('/Users/ertailuo/Research/CORA-2024/multiaff/data/supports/verts_supports_matlab.mat');
old_pts = pts;
load('/Users/ertailuo/Research/CORA-2024/multiaff/data/supports/verif_pts.mat');
assert (isequal(old_pts, pts), 'Not equal pts');

% Assert the directions are nomal vectors =================================
norm_sqr = dirs(:, 1).^2 + dirs(:, 2).^2;
ln = length(norm_sqr(:, 1));
for i = 1:ln
    assert(abs(norm_sqr(i) - 1) < 1e-7, "norms not normalized");
end

% Compute the f* values ===================================================
d1 = dirs(:, 1);
p1 = pts(:, 1);
prod1 = d1 .* p1;

d2 = dirs(:, 2);
p2 = pts(:, 2);
prod2 = d2 .* p2;

f_star = prod1 + prod2;

% Construct the diray name and result matfile storing both the support and
% the corresponding accuracy <--> min(|f1 - f*| / |f*|, |f2 - f*| / |f*|)
if startidx == 1
    diary_name = sprintf('%s%d%s', "multiaff/out/support/support_log_split_", ...
            split_num,'.txt');
    mat_name = sprintf('%s%d%s', "multiaff/out/support/support_split_", ...
            split_num, '_itv_dev.mat');
else
    diary_name = sprintf('%s%d%s%d%s', "multiaff/out/support/support_log_split_", ...
            split_num,'_start_idx_',startidx,'.txt');
    mat_name = sprintf('%s%d%s%d%s', "multiaff/out/support/support_split_", ...
            split_num, '_start_idx_',startidx, '_itv_dev.mat');
end

dir_num = ln;
result_intv = cell(ln, 1);
result_dev = cell(ln, 1);

diary(diary_name);
disp(['The optimal values of f* are: [', num2str(f_star'), ']']);

for i = startidx : dir_num
    dir = dirs(i, :)';
    disp(['Case for i = ', num2str(i)]);
    starttime = datetime('now');
    disp(['Start time of computing support for i = ', num2str(i), ...
        ' direction == [', num2str(dir'), '], with split num: ', num2str(split_num), ...
        ' is: ', char(starttime)]);

    val = supportFunc_(pZ_set, dir, 'range', 'split', split_num);

    endtime = datetime('now');
    disp(['End time of computing support for i = ', num2str(i), ...
        ' direction == [', num2str(dir'), '], with split num: ', num2str(split_num), ...
        ' is: ', char(endtime)]);
    
    % Store and display the result support interval
    result_intv{i} = val;
    inf = val.inf;
    sup = val.sup;
    disp(['Resulting support interval has inf: ', num2str(inf), ...
        ' sup: ', num2str(sup)]);

    % Compute the accuracy / deviation through 
    %   min(|f1 - f*| / |f*|, |f2 - f*| / |f*|)
    % The smaller, the closer
    temp_f_star = f_star(i);
    dev1 = abs(inf - temp_f_star) / abs(temp_f_star);
    dev2 = abs(sup - temp_f_star) / abs(temp_f_star);
    deviation = min(dev1, dev2);
    result_dev{i} = deviation;

    disp(['Resulting deviation is: ', num2str(deviation)]);
end

diary off;

% Dump the storeed support interval together with the deviation
save(mat_name, "result_intv", "result_dev");
clear;
end

