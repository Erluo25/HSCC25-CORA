function vardubin_exp(option)
% =========================================================================
%{
Shared params structure for easy copy past
-------------------------------------------
params.given_pZ = ;
params.pzfilename = ;
params.dirfilename = ;
params.origin_splits_num = ;
params.itct_split_num = ;
params.projdim = ;
params.option = ;
-------------------------------------------
params.rdpts_num = ;
params.figfname = ;
%}
% =========================================================================
if strcmp(option, "compare")
    load('/Users/ertailuo/Research/CORA-2024/multiaff/data/var_dubins/verification_vardubin_pz.mat');
    params.given_pZ = veri_pZ_set;
    params.pzfilename = "multiaff/data/var_dubins/struct_vardubins.mat";
    params.dirfilename = "multiaff/data/var_dubins/verts_var_dubins_matlab.mat";
    params.origin_splits_num = 20;
    params.itct_split_num = 20;
    params.projdim = [1 2];
    params.option = option;
    params.figfname = "../Projects/Multi-affine/Journal Figures/compare_60_ours_base.fig";
    mtaff_new_plot(params);
end

if strcmp(option, "test")
    load('/Users/ertailuo/Research/CORA-2024/multiaff/data/var_dubins/verification_vardubin_pz.mat');
    params.given_pZ = veri_pZ_set;
    params.pzfilename = "multiaff/data/var_dubins/struct_vardubins.mat";
    params.dirfilename = "multiaff/data/var_dubins/verts_var_dubins_matlab.mat";
    params.origin_splits_num = 20;
    params.itct_split_num = 20;
    params.projdim = [1 2];
    params.option = option;
    mtaff_new_plot(params);
end

if strcmp(option, "default")
    load('/Users/ertailuo/Research/CORA-2024/multiaff/data/var_dubins/verification_vardubin_pz.mat');
    params.given_pZ = veri_pZ_set;
    params.pzfilename = "multiaff/data/var_dubins/struct_vardubins.mat";
    params.dirfilename = "multiaff/data/var_dubins/verts_var_dubins_matlab.mat";
    params.origin_splits_num = 60;
    params.itct_split_num = -1;
    params.projdim = [1 2];
    params.option = option;
    params.rdpts_num = 100000;
    mtaff_new_plot(params);
end

if strcmp(option, "countor")
    load('/Users/ertailuo/Research/CORA-2024/multiaff/data/var_dubins/verification_vardubin_pz.mat');
    params.given_pZ = veri_pZ_set;
    params.pzfilename = "multiaff/data/var_dubins/struct_vardubins.mat";
    params.dirfilename = "multiaff/data/var_dubins/verts_var_dubins_matlab.mat";
    params.origin_splits_num = 20;
    params.itct_split_num = 20;
    params.projdim = [1 2];
    params.option = option;
    mtaff_new_plot(params);
end

end

