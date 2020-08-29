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

% dataset = 'warps_plane_trial11.mat';

% dataset = 'warps_cylinder1.mat';
% dataset = 'warps_cylinder2.mat';
% dataset = 'warps_cylinder3.mat';

show_plot = 0; % flag for showing the recovered shapes
decomp = 1;
method = struct;

method.method = 1; 
% 0 for local approach: locally select the shape parameters
% 1 for global approach: select the solution by maximizing the consistency
% 2 for least median: select the least median (not a two-view method)


%%% local approach
% method.measure = 'msa'; % minimum surface area
% method.measure = 'apap'; % as parallel as possible
method.measure = 'ln'; % least norm or least change of depth

method.visb = 1; % flag for applying visibility condition
method.both_view = 0; % flag for penalizing on both views


%%% global approach
method.solver = 'admm';
% method.solver = 'qp';
method.c1 = 1.0; % parallelism for both views
method.c2 = 0.1; % least depth change in both views
method.c3 = 1.0; % visibility for both views

% parameters for constructing the Laplacian
method.sigma = 0.2;
method.ratio = 1.0; % threshold = ratio * average squared distance


dataset = ['./warps_data/', dataset];
load(dataset);
pairs_num = length(qgth) - 1;

%3D Ground truth points and normalized image points
n = length(qgth);
for i=1:n
    q_n(2*(i-1)+1:2*(i-1)+2,:) = qgth{i}(1:2, :);
end
%visiblity matrix: remove points visible in less than 3 views
visb = ones(n, size(I1u, 2));

%%% GROUND TRUTH NORMALS %%%%%
% Ngth = create_gth_normals(Pgth,q_n,n);

%%% PARAMETERS %%%%%
par = 2e-3; % schwarzian parameter.. needs to be tuned (usually its something close to 1e-3)

% add noise
rng(2020);
pixel_noise = 5;
f = 500; num_p = length(I1u(1, :));
I1u(1, :) = randn(1, num_p) * pixel_noise / f + I1u(1, :);
I1v(1, :) = randn(1, num_p) * pixel_noise / f + I1v(1, :);
I2u = randn(pairs_num, num_p) * pixel_noise / f + I2u;
I2v = randn(pairs_num, num_p) * pixel_noise / f + I2v;


I1u = repmat(I1u(1, :), pairs_num, 1);
I1v = repmat(I1v(1, :), pairs_num, 1);


[I1u,I1v,I2u,I2v,J21a,J21b,J21c,J21d,J12a,J12b,J12c,J12d,H21uua,H21uub,H21uva,H21uvb,H21vva,H21vvb] = create_warps(I1u,I1v,I2u,I2v,visb,par);



view_id_list = 2:(1 + pairs_num);
% view_id_list = [2];

direct_substitute = 0;
second_view_id = 2;


err = 1e-20; % tolerance for imaginary part 
%%% (if a number's imaginary part is greater than err, then it's deemed as complex number)


num = 2; % 2 views

error_metric = zeros(4, length(view_id_list));
I1u = I1u(1, :);
I1v = I1v(1, :);

eq_coef = cell(1, pairs_num); % coefficients for the 6th degree polynomial

% coefficients for second cubic (f1 x1 + f2 = 0)
f1_coef = cell(1, pairs_num);
f2_coef = cell(1, pairs_num);

first_cubic_coef = cell(1, pairs_num); % coefficients for the first cubic
J21_all = cell(1, pairs_num); % collect a,b,c,d
t_all = cell(1, pairs_num); % collect t1, t2

tic
% collect coefficients
for view_id = 2:(1 + pairs_num)
    % keep only the first view and another view
    I2u_tem = I2u(view_id - 1, :); I2v_tem = I2v(view_id - 1, :);

    % coeff of k1       % coeff of k2        % constant term
    T1_12_k1 = -J21c(view_id - 1, :);   T1_12_k2 = -J21d(view_id - 1, :);    T1_12_c = (J12a(view_id - 1, :).*H21uva(view_id - 1, :) + J12c(view_id - 1, :).*H21uvb(view_id - 1, :));%(H21vvb./J21d)/2;
    T2_12_k1 = -J21a(view_id - 1, :);   T2_12_k2 = -J21b(view_id - 1, :);    T2_12_c = (J12b(view_id - 1, :).*H21uva(view_id - 1, :) + J12d(view_id - 1, :).*H21uvb(view_id - 1, :));%(H21uub./J21b)/2;

    % k1b = -T2_12 = a*k1 + b*k2 + t1;
    % k2b = -T1_12 = c*k1 + d*k2 + t2;
    a = -T2_12_k1; b = -T2_12_k2; c = -T1_12_k1; d = -T1_12_k2; t1 = -T2_12_c; t2 = -T1_12_c;
    
    J21_all{view_id - 1} = [a;b;c;d]; t_all{view_id - 1} = [t1; t2];
    
    e = 1+ I2u_tem.^2 + I2v_tem.^2; u = I2u_tem; v = I2v_tem;
    e1 = 1+ I1u.^2 + I1v.^2; u1 = I1u; v1 = I1v;
    
    [eq, f1, f2, r] = create_cubic_coefficients(a,b,c,d,t1,t2,e,e1,u,u1,v,v1);
    eq_coef{view_id - 1} = eq;
    f1_coef{view_id - 1} = f1;
    f2_coef{view_id - 1} = f2;
    first_cubic_coef{view_id - 1} = r;
end
toc

method.eq_coef = eq_coef;
method.f1_coef = f1_coef;
method.f2_coef = f2_coef;
method.first_cubic_coef = first_cubic_coef;
method.t_all = t_all;
method.J21_all = J21_all;


if direct_substitute

    % start solving all the cubics
    view_id = second_view_id;
    method.view_id = view_id;
    I2u_tem = I2u(view_id - 1, :); I2v_tem = I2v(view_id - 1, :);

    % solve x1 and x2
    eq = eq_coef{view_id - 1};
    f1 = f1_coef{view_id - 1};
    f2 = f2_coef{view_id - 1};
    [x1, x2] = solve_cubic(eq, f1, f2, err);
    %%% Note that we replace all complex numbers with NaN
    %%% there are 6 roots, so x1 and x2 have 6 rows of solution

    J12a_all = repmat(J12a(view_id - 1, :), 6, 1); J12b_all = repmat(J12b(view_id - 1, :), 6, 1);
    J12c_all = repmat(J12c(view_id - 1, :), 6, 1); J12d_all = repmat(J12d(view_id - 1, :), 6, 1);

    % on the first view
    k1 = J12a_all .* x1 + J12b_all .* x2;
    k2 = J12c_all .* x1 + J12d_all .* x2;

    % on the second view
    k1_ = x1 + repmat(t_all{view_id - 1}(1, :), 6, 1);
    k2_ = x2 + repmat(t_all{view_id - 1}(2, :), 6, 1);
    
    % compute the mask for selecting the desirable solutions
    mask = solution_selection(I1u, I1v, I2u_tem, I2v_tem, k1, k2, k1_, k2_, method);

    % gather k1 and k2 for both views
    k1 = k1(mask)'; % k1_ = k1_(mask)';
    k2 = k2(mask)'; % k2_ = k2_(mask)';
%     k1_all = [k1; k1_];
%     k2_all = [k2; k2_];
    count = 1;
    for view_id = view_id_list
        Ngth_tem = [Ngth(1:3, :); Ngth((view_id * 3 - 2):(view_id * 3), :)];
        Pgth_tem = [Pgth(1:3, :); Pgth((view_id * 3 - 2):(view_id * 3), :)];
        q_n = [qgth{1}(1:2, :); qgth{view_id}(1:2, :)];
        
        u_all = [I1u(1,:);I2u(view_id - 1, :)]; v_all = [I1v(1,:);I2v(view_id - 1, :)];

        k1_ = J12a(view_id - 1, :) .* k1 + J12b(view_id - 1, :) .* k2 + t_all{view_id - 1}(1, :);
        k2_ = J12c(view_id - 1, :) .* k1 + J12d(view_id - 1, :) .* k2 + t_all{view_id - 1}(2, :);

        k1_all = [k1; k1_];
        k2_all = [k2; k2_];
        
        % find normals on all surfaces N= [N1;N2;N3]
        N1 = k1_all; N2 = k2_all; N3 = 1-u_all.*k1_all-v_all.*k2_all;
        n = sqrt(N1.^2+N2.^2+N3.^2);
        N1 = N1./n ; N2 = N2./n; N3 = N3./n;

        N = [N1(:),N2(:),N3(:)]';
        N_res = reshape(N(:),3*num, length(u_all));

        % Integrate normals to find depth
        P_grid=calculate_depth(N_res,u_all,v_all,1e0);
        % compare with ground truth
        [P2,err_p] = compare_with_Pgth(P_grid,u_all,v_all,q_n,Pgth_tem);
        [N,err_n] = compare_with_Ngth(P2,q_n,Ngth_tem);
        error_metric(1:2, count) = mean(err_p');
        error_metric(3:4, count) = mean(err_n');

        if show_plot
            figure();
    %         ax = gca;
    %         ax.ZAxis.TickLabelFormat = '%.2f';        

            draw_surface(Pgth_tem(1:3, :), 'g')

            hold on
            draw_surface(P2(1:3, :), 'r')
            hold off
            figure();
            draw_surface(Pgth_tem(4:6, :), 'g')
            hold on
            draw_surface(P2(4:6, :), 'r')
            hold off
        end
        count = count + 1;
    end
    %     disp('Frobenius norm of errors of the computed normals')
    %     disp(norm(-N_res - Ngth_tem))
else
    
    T_poly = zeros(1, length(view_id_list));
    T_sel = zeros(1, length(view_id_list));
    T_norm = zeros(1, length(view_id_list));
    
    count = 1;
    % start solving all the cubics
    for view_id = view_id_list
        method.view_id = view_id;
        I2u_tem = I2u(view_id - 1, :); I2v_tem = I2v(view_id - 1, :);
        t1 = t_all{view_id - 1}(1, :);
        t2 = t_all{view_id - 1}(2, :);

        Ngth_tem = [Ngth(1:3, :); Ngth((view_id * 3 - 2):(view_id * 3), :)];
        Pgth_tem = [Pgth(1:3, :); Pgth((view_id * 3 - 2):(view_id * 3), :)];
        q_n = [qgth{1}(1:2, :); qgth{view_id}(1:2, :)];

        tic
        if decomp
            num_p = length(t1);
            T = eye(3);
            k1 = zeros(2, num_p);
            k2 = zeros(2, num_p);
            k1_ = zeros(2, num_p);
            k2_ = zeros(2, num_p);
            
            for i =1:num_p
                T_x = T; T_y = T;
                x = [I1u(i); I1v(i)]; y = [I2u_tem(i); I2v_tem(i)];
                T_x(3, 1:2) = x';
                T_y(3, 1:2) = -y';
                T_k = [J21a(view_id - 1, i), J21b(view_id - 1, i), t1(i);
                        J21c(view_id - 1, i), J21d(view_id - 1, i), t2(i);
                        0, 0, 1];
                W = T_y * T_k * T_x;
                
                [n1_, n2_] = compute_normal(W');
                k_second = [n1_(1:2) / ([y; 1]' * n1_), n2_(1:2) / ([y; 1]' * n2_)];
                k1_(:, i) = k_second(1, :)';
                k2_(:, i) = k_second(2, :)';
                
                n1 = W \ n1_;
                n2 = W \ n2_;
                k_first = [n1(1:2) / ([x; 1]' * n1), n2(1:2) / ([x; 1]' * n2)];
                k1(:, i) = k_first(1, :)';
                k2(:, i) = k_first(2, :)';
            end
            
            tem = nan(4, num_p);
            k1 = [tem; k1];
            k2 = [tem; k2];
            k1_ = [tem; k1_];
            k2_ = [tem; k2_];
            
            
        else
            % solve x1 and x2
            eq = eq_coef{view_id - 1};
            f1 = f1_coef{view_id - 1};
            f2 = f2_coef{view_id - 1};
            
            [x1, x2] = solve_cubic(eq, f1, f2, err);
            
            %%% Note that we replace all complex numbers with NaN
            %%% there are 6 roots, so x1 and x2 have 6 rows of solution

            J12a_all = repmat(J12a(view_id - 1, :), 6, 1); J12b_all = repmat(J12b(view_id - 1, :), 6, 1);
            J12c_all = repmat(J12c(view_id - 1, :), 6, 1); J12d_all = repmat(J12d(view_id - 1, :), 6, 1);

            % on the first view
            k1 = J12a_all .* x1 + J12b_all .* x2;
            k2 = J12c_all .* x1 + J12d_all .* x2;

            % on the second view
            k1_ = x1 + repmat(t1, 6, 1);
            k2_ = x2 + repmat(t2, 6, 1);
        end
        T_poly(count) = toc;
        
        % compute the mask for selecting the desirable solutions
        tic
        mask = solution_selection(I1u, I1v, I2u_tem, I2v_tem, k1, k2, k1_, k2_, method);
        T_sel(count) = toc;
        
        % gather k1 and k2 for both views
        k1_all = [k1(mask)'; k1_(mask)'];
        k2_all = [k2(mask)'; k2_(mask)'];

        u_all = [I1u(1,:);I2u_tem]; v_all = [I1v(1,:);I2v_tem];

        % find normals on all surfaces N= [N1;N2;N3]
        N1 = k1_all; N2 = k2_all; N3 = 1-u_all.*k1_all-v_all.*k2_all;
        n = sqrt(N1.^2+N2.^2+N3.^2);
        N1 = N1./n ; N2 = N2./n; N3 = N3./n;

        N = [N1(:),N2(:),N3(:)]';
        N_res = reshape(N(:),3*num, size(u_all, 2));

        tic
        % Integrate normals to find depth
        P_grid = calculate_depth(N_res,u_all,v_all,1e0);
        
        % compare with ground truth
        [P2,err_p] = compare_with_Pgth(P_grid,u_all,v_all,q_n,Pgth_tem);
        [N,err_n] = compare_with_Ngth(P2,q_n,Ngth_tem);
        T_norm(count) = toc;
        
        error_metric(1:2, count) = mean(err_p');
        error_metric(3:4, count) = mean(err_n');

        if show_plot
            figure();
    %         ax = gca;
    %         ax.ZAxis.TickLabelFormat = '%.2f';        

            draw_surface(Pgth_tem(1:3, :), 'g')

            hold on
            draw_surface(P2(1:3, :), 'r')
            hold off
            figure();
            draw_surface(Pgth_tem(4:6, :), 'g')
            hold on
            draw_surface(P2(4:6, :), 'r')
            hold off
        end
    %     disp('Frobenius norm of errors of the computed normals')
%         disp(norm(-N_res - Ngth_tem))
        count = count + 1;
    end
end



if length(view_id_list) > 1
    disp('Errors for the points')
    mean(error_metric(1:2, :)')
    disp('Errors for the normals (first and second)')
    mean(error_metric(3:4, :)')
    
    disp('Errors for the normals (all)')
    mean([mean(error_metric(3, :)), error_metric(4, :)])
    
else
    disp('Errors for the points')
    error_metric(1:2, :)'
    disp('Errors for the normals')
    error_metric(3:4, :)'
end

