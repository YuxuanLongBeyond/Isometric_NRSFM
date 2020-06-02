% example script
clear all;close all;

% add libraries
addpath(genpath('BBS_NOOMP'));
addpath(genpath('tbxmanager'));
tbxmanager restorepath

addpath('SfTv0_3');
addpath(genpath('gloptipoly3'));
addpath(genpath('SeDuMi_1_3'));
addpath(genpath('schwarps'));
addpath(genpath('sparseinv'));
addpath(genpath('utils'));
% dataset = 'Kinect_paper.mat';
% dataset = 'warps_tshirt.mat';
% dataset = 'warps_plane1.mat';
% dataset = 'warps_plane2.mat';
dataset = 'warps_cylinder1.mat';
% dataset = 'warps_cylinder2.mat';
% dataset = 'warps_cylinder3.mat';


dataset = ['./warps_data/', dataset];
load(dataset);

%3D Ground truth points and normalized image points
num = length(qgth);
for i=1:num
    q_n(2*(i-1)+1:2*(i-1)+2,:) = qgth{i}(1:2, :);
end
%visiblity matrix: remove points visible in less than 3 views
visb = ones(num, size(Pgth, 2));

%%% GROUND TRUTH NORMALS %%%%%
Ngth = create_gth_normals(Pgth,q_n,num);

%%% PARAMETERS %%%%%
par = 2e-3; % schwarzian parameter.. needs to be tuned (usually its something close to 1e-3)
% grid = 1;
% 
% %%%%% GRID OF POINTS %%%%
% if grid %max(sum(visb)) == num
%     % make a grid
%     [I1u,I1v,I2u,I2v,visb] = create_grid(q_n,visb,20);

% else
%     % point-wise
%     I1u = repmat(q_n(1,:),num-1,1);
%     I1v = repmat(q_n(2,:),num-1,1);
%     I2u = q_n(3:2:2*num,:);
%     I2v = q_n(4:2:2*num,:);
% end
%%%%% SCHWARZIAN WARPS %%%%%
% [I1u,I1v,I2u,I2v,J21a,J21b,J21c,J21d,J12a,J12b,J12c,J12d,H21uua,H21uub,H21uva,H21uvb,H21vva,H21vvb] = create_warps(I1u,I1v,I2u,I2v,visb,par);
% create schwarzian warps for the dataset

% Christoffel Symbols (see equation 15 in the paper)
%T1 = [-2*k1 -k2;-k2 0];
%T2 = [0 -k1;-k1 -2*k2];
% Christoffel Symbols change of variable  (see equation 10 in the paper) written in terms of k1b and k2b

if size(I1u, 1) == 1
    I1u = repmat(I1u, num-1,1);
    I1v = repmat(I1v ,num-1,1);
end

% coeff of k1       % coeff of k2        % constant term
T1_12_k1 = -J21c;   T1_12_k2 = -J21d;    T1_12_c = (J12a.*H21uva + J12c.*H21uvb);%(H21vvb./J21d)/2;
T2_12_k1 = -J21a;   T2_12_k2 = -J21b;    T2_12_c = (J12b.*H21uva + J12d.*H21uvb);%(H21uub./J21b)/2;

% k1b = -T2_12 = a*k1 + b*k2 + t1;
% k2b = -T1_12 = c*k1 + d*k2 + t2;
a = -T2_12_k1; b = -T2_12_k2; c = -T1_12_k1; d = -T1_12_k2; t1 = -T2_12_c; t2 = -T1_12_c;

e = 1+ I2u.^2 + I2v.^2; u = I2u; v = I2v;
e1 = 1+ I1u.^2 + I1v.^2; u1 = I1u; v1 = I1v;

% create sum of squares polynomial
eq = create_polynomial_coefficients(a,b,c,d,t1,t2,e,e1,u,u1,v,v1);
% minimise it to obtain depth derivatives
tic
res = solve_polynomial(eq);
toc
% recover first order derivatives on rest of the surfaces
k1_all = [res(:,1)';a.*repmat(res(:,1)',num-1,1) + b.*repmat(res(:,2)',num-1,1) + t1];
k2_all = [res(:,2)';c.*repmat(res(:,1)',num-1,1) + d.*repmat(res(:,2)',num-1,1) + t2];

idx = find(visb(1,:)==0);
for i = 1: length(idx)
    id = find(visb(1:end,idx(i))>0);
    I2u(id(1)-1,idx(i)) = I1u(1,idx(i)); I2v(id(1)-1,idx(i)) = I1v(1,idx(i)); I1u(:,idx(i)) = 0; I1v(:,idx(i)) = 0;
    k1_all(id(1),idx(i)) = k1_all(1,idx(i)); k2_all(id(1),idx(i)) = k2_all(1,idx(i)); k1_all(1,idx(i)) = 0; k2_all(1,idx(i)) = 0;
end

u_all = [I1u(1,:);I2u]; v_all = [I1v(1,:);I2v];

% find normals on all surfaces N= [N1;N2;N3]
N1 = k1_all; N2 = k2_all; N3 = 1-u_all.*k1_all-v_all.*k2_all;
n = sqrt(N1.^2+N2.^2+N3.^2);
N1 = N1./n ; N2 = N2./n; N3 = N3./n;

N = [N1(:),N2(:),N3(:)]';
N_res = reshape(N(:),3*num,length(u_all));

% find indices with no solution
idx = find(res(:,1)==0);
N_res(:,idx) = []; u_all(:,idx) = []; v_all(:,idx) = [];

% Integrate normals to find depth
P_grid=calculate_depth(N_res,u_all,v_all,1e0);

% compare with ground truth
[P2,err_p] = compare_with_Pgth(P_grid,u_all,v_all,q_n,Pgth);
[N,err_n] = compare_with_Ngth(P2,q_n,Ngth);

% % plot results
% for i=1:size(u_all,1)
%      figure(i)
%     plot3(Pgth(3*(i-1)+1,:),Pgth(3*(i-1)+2,:),Pgth(3*(i-1)+3,:),'go');
%     hold on;
%       plot3(P2(3*(i-1)+1,:),P2(3*(i-1)+2,:),P2(3*(i-1)+3,:),'ro');
%       %quiver3(P2(3*(i-1)+1,:),P2(3*(i-1)+2,:),P2(3*(i-1)+3,:),N(3*(i-1)+1,:),N(3*(i-1)+2,:),N(3*(i-1)+3,:));
%       %quiver3(Pgth(3*(i-1)+1,:),Pgth(3*(i-1)+2,:),Pgth(3*(i-1)+3,:),Ngth(3*(i-1)+1,:),Ngth(3*(i-1)+2,:),Ngth(3*(i-1)+3,:));
%     hold off;
%     axis equal;
% end
mean(err_p');
% mean(err_n');
tem = mean(err_n');
disp('Average normal error on first frame')
mean(tem(1))
disp('Average normal error on other frames')
mean(tem(2:end))
