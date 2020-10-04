clear all
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

data = {'cat', 'rug_trun', 'Kinect_paper', 'tshirt', 'warps_cylinder3'};
dataset = data{5};
load(['./warps_data/', dataset, '.mat']);
pairs_num = length(qgth) - 1;
num = length(qgth);
for i=1:num
    q_n(2*(i-1)+1:2*(i-1)+2,:) = qgth{i}(1:2, :);
end
visb = ones(num, size(q_n, 2));
Ngth = create_gth_normals(Pgth,q_n,num, visb);
% [I1u,I1v,I2u,I2v,visb] = create_grid(q_n,visb,20);
I1u = repmat(q_n(1,:),num-1,1);
I1v = repmat(q_n(2,:),num-1,1);
I2u = q_n(3:2:2*num,:);
I2v = q_n(4:2:2*num,:);

par = 2e-3;
[I1u,I1v,I2u,I2v,J21a,J21b,J21c,J21d,J12a,J12b,J12c,J12d,H21uua,H21uub,H21uva,H21uvb,H21vva,H21vvb] = create_warps(I1u,I1v,I2u,I2v,visb,par, 1);


save([dataset, '_nongrid.mat'], 'H21uua', 'H21uub', 'H21uva', 'H21uvb', 'H21vva', 'H21vvb', ...
    'I1u', 'I1v', 'I2u', 'I2v', 'J21a', 'J21b', 'J21c', 'J21d', 'J12a', 'J12b', 'J12c', 'J12d', ...
    'Ngth', 'Pgth', 'qgth');
