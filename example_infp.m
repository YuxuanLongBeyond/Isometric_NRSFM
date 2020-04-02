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

%%%% INPUTS %%%%%
% load tshirt.mat

%3D Ground truth points and normalized image points
% num = length(scene.m);
% for i=1:num
%     Pgth(3*(i-1)+1:3*(i-1)+3,:) = scene.Pgth(i).P; 
%     q_n(2*(i-1)+1:2*(i-1)+2,:) = scene.m(i).m(1:2,:);
% end
%visiblity matrix: remove points visible in less than 3 views
% visb = ones(num,length(scene.Pgth(1).P(1,:)));

%%% GROUND TRUTH NORMALS %%%%%
% Ngth = create_gth_normals(Pgth,q_n,num);

% %%% PARAMETERS %%%%%
% par = 2e-3; % schwarzian parameter.. needs to be tuned (usually its something close to 1e-3)
% grid = 0;

%%%%% GRID OF POINTS %%%%
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

% load('Kinect_paper.mat')
% load('warps_tshirt.mat')
% load('warps_plane_vertical.mat')
% load('warps_plane_horizontal.mat')
load('warps_cylinder.mat')

% select one view
view_id = 3; % from 1 to 10

% solution selection by methods:
method = struct;
% method.method = 0; % select the solution by seeking the minimum absolute value
method.method = 1; % select the solution by exploring the graph Laplacian
% method.method = 3;
method.sigma = 1;
method.ratio = 1; % threshold = ratio * average squared distance
method.solver = 'qp1'; % no inequality constraint
% method.solver = 'qp2'; % with inequality constraint
% method.solver = 'admm'; % with l1 penalty (to replace inequality constraint)

% method.solver = 'irqp'; 
num = 2;

% keep only the first view and another view
H21uua = H21uua(view_id - 1, :);
H21uub = H21uub(view_id - 1, :);
H21uva = H21uva(view_id - 1, :);
H21uvb = H21uvb(view_id - 1, :);
H21vva = H21vva(view_id - 1, :);
H21vvb = H21vvb(view_id - 1, :);
I1u = I1u(1, :);
I1v = I1v(1, :);
I2u = I2u(view_id - 1, :);
I2v = I2v(view_id - 1, :);
J12a = J12a(view_id - 1, :);
J12b = J12b(view_id - 1, :);
J12c = J12c(view_id - 1, :);
J12d = J12d(view_id - 1, :);
J21a = J21a(view_id - 1, :);
J21b = J21b(view_id - 1, :);
J21c = J21c(view_id - 1, :);
J21d = J21d(view_id - 1, :);
Ngth = [Ngth(1:3, :); Ngth((view_id * 3 - 2):(view_id * 3), :)];
Pgth = [Pgth(1:3, :); Pgth((view_id * 3 - 2):(view_id * 3), :)];
q_n = [qgth{1}(1:2, :); qgth{view_id}(1:2, :)];

% Christoffel Symbols (see equation 15 in the paper)
%T1 = [-2*k1 -k2;-k2 0];
%T2 = [0 -k1;-k1 -2*k2];
% Christoffel Symbols change of variable  (see equation 10 in the paper) written in terms of k1b and k2b

% coeff of k1       % coeff of k2        % constant term
T1_12_k1 = -J21c;   T1_12_k2 = -J21d;    T1_12_c = (J12a.*H21uva + J12c.*H21uvb);%(H21vvb./J21d)/2;
T2_12_k1 = -J21a;   T2_12_k2 = -J21b;    T2_12_c = (J12b.*H21uva + J12d.*H21uvb);%(H21uub./J21b)/2;

% k1b = -T2_12 = a*k1 + b*k2 + t1;
% k2b = -T1_12 = c*k1 + d*k2 + t2;
a = -T2_12_k1; b = -T2_12_k2; c = -T1_12_k1; d = -T1_12_k2; t1 = -T2_12_c; t2 = -T1_12_c;

e = 1+ I2u.^2 + I2v.^2; u = I2u; v = I2v;
e1 = 1+ I1u.^2 + I1v.^2; u1 = I1u; v1 = I1v;

[eq, f1, f2] = create_cubic_coefficients(a,b,c,d,t1,t2,e,e1,u,u1,v,v1);
% minimise it to obtain depth derivatives
err = 1e-10;
[x1, x2] = solve_cubic(eq, f1, f2, err);


J12a_all = repmat(J12a, 6, 1); J12b_all = repmat(J12b, 6, 1);
J12c_all = repmat(J12c, 6, 1); J12d_all = repmat(J12d, 6, 1);
I1u_all = repmat(I1u, 6, 1); I1v_all = repmat(I1v, 6, 1);

% on the first view
k1 = J12a_all .* x1 + J12b_all .* x2;
k2 = J12c_all .* x1 + J12d_all .* x2;

% on the second view
k1_ = x1 + repmat(t1, 6, 1);
k2_ = x2 + repmat(t2, 6, 1);

% n3 = 1 - k1 .* I1u_all - k2 .* I1v_all;
% n3_ = 1 - k1_ .* repmat(I2u, 6, 1) - k2_ .* repmat(I2v, 6, 1);
% mask = (n3 < 0) | (n3_ < 0); % the normal should point outwards under visible condition
% k1(mask) = NaN; k2(mask) = NaN; 
% x1(mask) = NaN; x2(mask) = NaN; 
mask = solution_selection(I1u, I1v, I2u, I2v, k1, k2, k1_, k2_, method);
k1_all = [k1(mask)'; k1_(mask)'];
k2_all = [k2(mask)'; k2_(mask)'];
% idx = find(visb(1,:)==0);
% for i = 1: length(idx)
%     id = find(visb(1:end,idx(i))>0);
%     I2u(id(1)-1,idx(i)) = I1u(1,idx(i)); I2v(id(1)-1,idx(i)) = I1v(1,idx(i)); I1u(:,idx(i)) = 0; I1v(:,idx(i)) = 0;
%     k1_all(id(1),idx(i)) = k1_all(1,idx(i)); k2_all(id(1),idx(i)) = k2_all(1,idx(i)); k1_all(1,idx(i)) = 0; k2_all(1,idx(i)) = 0;
% end

u_all = [I1u(1,:);I2u]; v_all = [I1v(1,:);I2v];

% find normals on all surfaces N= [N1;N2;N3]
N1 = k1_all; N2 = k2_all; N3 = 1-u_all.*k1_all-v_all.*k2_all;
n = sqrt(N1.^2+N2.^2+N3.^2);
N1 = N1./n ; N2 = N2./n; N3 = N3./n;

N = [N1(:),N2(:),N3(:)]';
N_res = reshape(N(:),3*num,length(u_all));

% find indices with no solution
% idx = find(res(:,1)==0);
% N_res(:,idx) = []; u_all(:,idx) = []; v_all(:,idx) = [];

% Integrate normals to find depth
P_grid=calculate_depth(N_res,u_all,v_all,1e0);

% compare with ground truth
[P2,err_p] = compare_with_Pgth(P_grid,u_all,v_all,q_n,Pgth);
[N,err_n] = compare_with_Ngth(P2,q_n,Ngth);

% plot results
for i=1:size(u_all,1)
     figure(i)
    plot3(Pgth(3*(i-1)+1,:),Pgth(3*(i-1)+2,:),Pgth(3*(i-1)+3,:),'go');
    hold on;
      plot3(P2(3*(i-1)+1,:),P2(3*(i-1)+2,:),P2(3*(i-1)+3,:),'ro');
      %quiver3(P2(3*(i-1)+1,:),P2(3*(i-1)+2,:),P2(3*(i-1)+3,:),N(3*(i-1)+1,:),N(3*(i-1)+2,:),N(3*(i-1)+3,:));
      %quiver3(Pgth(3*(i-1)+1,:),Pgth(3*(i-1)+2,:),Pgth(3*(i-1)+3,:),Ngth(3*(i-1)+1,:),Ngth(3*(i-1)+2,:),Ngth(3*(i-1)+3,:));
    hold off;
%     axis equal;
end
mean(err_p')
mean(err_n')
