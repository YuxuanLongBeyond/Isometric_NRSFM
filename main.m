% example script
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

use_warp = 1; % if warp is already conntained in the data
% schwarp = 0;

pixel_noise = 0;
f = 500;  % default focal length

show_plot = 0; % flag for showing the recovered shapes

grid = 1;
if strcmp(dataset(1:5), 'warps')
    grid = 0;
end


%% test two view methods
decomp = 1;
method = struct;

choice = 0; 
% 0 for local approach: locally select the shape parameters
% 1 for global approach: select the solution by maximizing the consistency
% 2 for least median: select the least median (not a two-view method)
% 3 for least dot product

%%% local approach
% measure = 'msa'; % minimum surface area
% measure = 'apap'; % as parallel as possible
measure = 'ln'; % least norm or least change of depth

use_visb = 1; % flag for applying visibility condition

%%% global approach
solver = 'admm';
% solver = 'qp';

[error_metric, T_poly, T_sel, T_norm] = test_two_view(dataset, pixel_noise, f, show_plot, decomp, choice, measure, use_visb, solver, grid, use_warp);
disp('Average time taken to reconstruct one pair of views:')
mean(T_poly) + mean(T_sel) + mean(T_norm)
err_n_ln = error_metric(3, :);
err_d_ln = error_metric(1, :);
save([dataset(1:end-4), '_two_view_ln_schwarp.mat' ], 'err_n_ln', 'err_d_ln')

%% test multiple view methods
% solver = 'infP';
% solver = 'iso';
% solver = 'polyH';
solver = 'fastDiffH';

views_num = 10; % should be greater than 2
% views_num = 'all';

tic
[err_n, err_p] = test_multiple_view(dataset, pixel_noise, f, solver, grid, show_plot, use_warp, views_num);
toc

%% test mode for multiple views

% solver = 'infP';
% solver = 'iso';
% solver = 'polyH';
solver = 'fastDiffH';

for views_num = 3:10
    [err_n, err_p] = test_multiple_view(dataset, pixel_noise, f, solver, grid, show_plot, use_warp, views_num);
    err_n1(views_num) = mean(err_n(1, :));
    err_d1(views_num) = mean(err_p(1, :));
    err_n_all(views_num) = mean(mean(err_n'));
    err_d_all(views_num) = mean(mean(err_p'));
end
save([dataset(1:end-4), '_multiple_view_', solver, '.mat'], 'err_n1', 'err_d1', 'err_n_all', 'err_d_all')

