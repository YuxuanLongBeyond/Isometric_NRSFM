clear all;close all;

% add libraries
addpath(genpath('BBS'));
addpath(genpath('tbxmanager'));
tbxmanager restorepath

addpath('SfTv0_3');
addpath(genpath('gloptipoly3'));
addpath(genpath('SeDuMi_1_3'));
addpath(genpath('schwarps'));
addpath(genpath('sparseinv'));
addpath(genpath('utils'));
addpath(genpath('l1magic'));


% dataset = 'Kinect_paper_nongrid.mat';
% dataset = 'rug_trun_nongrid.mat';
% dataset = 'cat_nongrid.mat';
% dataset = 'tshirt_nongrid.mat';

% dataset = 'Kinect_paper.mat';
% dataset = 'rug_trun.mat';
% dataset = 'cat.mat';
dataset = 'tshirt.mat';

% dataset = 'warps_plane1.mat';
% dataset = 'warps_plane2.mat';

% dataset = 'warps_plane_trial11.mat';

% dataset = 'warps_cylinder1.mat';
% dataset = 'warps_cylinder2.mat';
% dataset = 'warps_cylinder3.mat';


choice = 0; 
% 0 for local approach: locally select the shape parameters
% 1 for global approach: select the solution by maximizing the consistency
% 2 for least median: select the least median (not a two-view method)
% 3 for least dot product

%%% local approach
% measure = 'msa'; % minimum surface area
% measure = 'apap'; % as parallel as possible
measure = 'ln'; % least norm or least change of depth

%%% global approach
% solver = 'admm';
solver = 'qp';

use_gth = 1;

use_warp = 0;
degen_filter = 0;

pixel_noise = 0;
grid = 1;
grid_size = 20;

show_plot = 0;
show_im = 0;

frame1 = 1;
load(['./warps_data/', dataset], 'qgth'); frame_num = length(qgth);

for frame2 = 2:frame_num
    [error_map1, error_map2, err_n, err_p] = two_view_nrsfm(dataset, frame1, frame2, pixel_noise, choice, measure, solver, grid, grid_size, use_warp, degen_filter, use_gth, show_plot, show_im);
    error_n1_raw(frame2 - 1) = mean(error_map1); error_n2_raw(frame2 - 1) = mean(error_map2);
    error_n1(frame2 - 1) = mean(err_n(1, :)); error_n2(frame2 - 1) = mean(err_n(2, :));
    error_p1(frame2 - 1) = mean(err_p(1, :)); error_p2(frame2 - 1) = mean(err_p(2, :));
end
disp('Raw average shape error of first frame')
mean(error_n1_raw)
disp('Raw average shape error of second frame')
mean(error_n2_raw)
disp('Average shape error of first frame')
mean(error_n1)
disp('Average shape error of second frame')
mean(error_n2)

disp('Average depth error of first frame')
mean(error_p1)
disp('Average depth error of second frame')
mean(error_p2)