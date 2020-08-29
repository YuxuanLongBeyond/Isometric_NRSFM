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
% dataset = 'Kinect_paper.mat';
dataset = 'warps_tshirt.mat';
% dataset = 'warps_plane1.mat';
% dataset = 'warps_plane2.mat';
% dataset = 'warps_cylinder1.mat';
% dataset = 'warps_cylinder2.mat';
% dataset = 'warps_cylinder3.mat';

fast_solve = 1;

dataset = ['./warps_data/', dataset];
load(dataset);
pairs_num = length(qgth) - 1;
%3D Ground truth points and normalized image points
num = length(qgth);
for i=1:num
    q_n(2*(i-1)+1:2*(i-1)+2,:) = qgth{i}(1:2, :);
end
%visiblity matrix: remove points visible in less than 3 views
visb = ones(num, size(I1u, 2));

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

% add noise
rng(2020);
pixel_noise = 0;
f = 500; num_p = length(I1u(1, :));
I1u(1, :) = randn(1, num_p) * pixel_noise / f + I1u(1, :);
I1v(1, :) = randn(1, num_p) * pixel_noise / f + I1v(1, :);
I2u = randn(num - 1, num_p) * pixel_noise / f + I2u;
I2v = randn(num - 1, num_p) * pixel_noise / f + I2v;


%%%%% SCHWARZIAN WARPS %%%%%
I1u = repmat(I1u(1, :), pairs_num, 1);
I1v = repmat(I1v(1, :), pairs_num, 1);
[I1u,I1v,I2u,I2v,J21a,J21b,J21c,J21d,J12a,J12b,J12c,J12d,H21uua,H21uub,H21uva,H21uvb,H21vva,H21vvb] = create_warps(I1u,I1v,I2u,I2v,visb,par);
% [~,~,~,I1u,I1v,I2u,I2v,J21a,J21b,J21c,J21d,J12a,J12b,J12c,J12d,H21uua,H21uub,H21uva,H21uvb,H21vva,H21vvb] = create_tshirt_dataset(idx,length(idx), scene, par);
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

tic

num_p = length(a(1, :));
N_res = zeros(3 * num, num_p);
T = eye(3);
if fast_solve
    for i = 1:num_p

        N = zeros(pairs_num, 3);
        W_set = cell(1, pairs_num);
        n1_all = zeros(3, pairs_num);
        n2_all = zeros(3, pairs_num);
        
        for j = 1:pairs_num

            T_x = T; T_y = T;
            x = [I1u(j, i); I1v(j, i)]; y = [I2u(j, i); I2v(j, i)];
            T_x(3, 1:2) = x';
            T_y(3, 1:2) = -y';
            T_k = [J21a(j, i), J21b(j, i), t1(j, i);
                    J21c(j, i), J21d(j, i), t2(j, i);
                    0, 0, 1];
            W = T_y * T_k * T_x;
            W_set{j} = W;

%             [R, t, n, d] = decompose(W',y,x);
            [n1, n2] = compute_normal(inv(W'));

            N(j, :) = cross(n1, n2);
            n1_all(:, j) = n1;
            n2_all(:, j) = n2;
        end
        [U,S,V] = svd(N);
        n = V(:, end);
        N_res(1:3, i) = n;  
        
        iter_max = 1;
        col_n = zeros(pairs_num, 3);
        for z = 1:iter_max
            P = zeros(3 * pairs_num, 3);
            for k = 1:pairs_num
                n1 = n1_all(:, k);
                n2 = n2_all(:, k);

                if abs(n' * n1) > abs(n' * n2)
                    P((3 * (k - 1) + 1):(3 * k), :) = eye(3) - n1 * n1';
                    col_n(k, :) = n1;
                else
                    P((3 * (k - 1) + 1):(3 * k), :) = eye(3) - n2 * n2';
                    col_n(k, :) = n2;
                end
            end

%             [U,S,V] = svd(P);
%             n = V(:, end);

            iter_max2 = 10;
            w = 1;
            for k = 1:iter_max2
                [U,S,V] = svd(sqrt(w) .* P);
                n = V(:, end);

                tem = P * n;
                w = 1 ./ sqrt(tem(1:3:end) .^ 2 + tem(2:3:end) .^ 2 + tem(3:3:end) .^ 2);
                w = [w, w, w]';
                w = w(:);
            end
        end
%         n = mean(col_n)'; n = n / norm(n);
        
        N_res(1:3, i) = n;  

        for j = 2:num
            n_ = W_set{j - 1} * n;
            N_res((3 * (j - 1) + 1):(3 * j), i) = n_ / norm(n_);
        end
    end
else
    for i = 1:num_p
        W_set = cell(1, pairs_num);

        coef_40 = 0; coef_31 = 0; coef_22 = 0;
        coef_13 = 0; coef_04 = 0; coef_30 = 0;
        coef_21 = 0; coef_12 = 0; coef_03 = 0;
        coef_20 = 0; coef_11 = 0; coef_02 = 0;
        coef_10 = 0; coef_01 = 0; coef_00 = 0;
        for j = 1:pairs_num

            T_x = T; T_y = T;
            x = [I1u(j, i); I1v(j, i)]; y = [I2u(j, i); I2v(j, i)];
            T_x(3, 1:2) = x';
            T_y(3, 1:2) = -y';
            T_k = [J21a(j, i), J21b(j, i), t1(j, i);
                    J21c(j, i), J21d(j, i), t2(j, i);
                    0, 0, 1];
            W = T_y * T_k * T_x;
            W_set{j} = W;
            H = (W')^(-1);
            
            sigma = svd(H);
            H = H / sigma(2);
            
            S = H' * H - eye(3);
            
            coef_40 = coef_40 + S(2, 2) ^ 2 + S(3, 3) ^ 2;
            coef_31 = coef_31 - 4 * S(1, 2) * S(2, 2);
            coef_22 = coef_22 + 2 * S(1, 1) * S(2, 2) + 4 * S(1, 2) ^ 2;
            coef_13 = coef_13 - 4 * S(1, 1) * S(1, 2);
            coef_04 = coef_04 + S(1, 1) ^ 2 + S(3, 3) ^ 2;
            coef_30 = coef_30 - 4 * S(1, 3) * S(3, 3);
            coef_03 = coef_03 - 4 * S(2, 3) * S(3, 3);
            coef_20 = coef_20 + 2 * S(1, 1) * S(3, 3) + 4 * S(1, 3) ^ 2;
            coef_02 = coef_02 + 2 * S(2, 2) * S(3, 3) + 4 * S(2, 3) ^ 2;
            coef_10 = coef_10 - 4 * S(1, 1) * S(1, 3);
            coef_01 = coef_01 - 4 * S(2, 2) * S(2, 3);
            coef_00 = coef_00 + S(1, 1) ^ 2 + S(2, 2) ^ 2;
        end
        eq = [coef_40; coef_31; coef_22; coef_13; coef_04; coef_30; coef_03; coef_20; coef_02; coef_10; coef_01; coef_00];
        n = [solve_normal(eq), 1]';
        n = n / norm(n);
        N_res(1:3, i) = n;
        
        for j = 2:num
            n_ = W_set{j - 1} * n;
            N_res((3 * (j - 1) + 1):(3 * j), i) = n_ / norm(n_);
        end
    end
end
toc

u_all = [I1u(1,:);I2u]; v_all = [I1v(1,:);I2v];

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
disp('Average normal error on all frames')
mean(tem)
