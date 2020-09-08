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
decomp = 0;
method = struct;

choice = 0; 
% 0 for local approach: locally select the shape parameters
% 1 for global approach: select the solution by maximizing the consistency
% 2 for least median: select the least median (not a two-view method)


%%% local approach
% measure = 'msa'; % minimum surface area
% measure = 'apap'; % as parallel as possible
measure = 'ln'; % least norm or least change of depth

use_visb = 1; % flag for applying visibility condition

%%% global approach
% solver = 'admm';
solver = 'qp';

direct_substitute = 0;
second_view_id = 2;
[error_metric, T_poly, T_sel, T_norm] = test_two_view(dataset, pixel_noise, f, show_plot, decomp, choice, measure, use_visb, solver, direct_substitute, second_view_id, grid, use_warp);
disp('Average time taken to reconstruct one pair of views:')
mean(T_poly) + mean(T_sel) + mean(T_norm)


%% test multiple view methods
% solver = 'infP';
% solver = 'iso';
% solver = 'polyH';
solver = 'fastDiffH';

tic
[err_n, err_p] = test_multiple_view(dataset, pixel_noise, f, solver, grid, show_plot, use_warp);
toc
