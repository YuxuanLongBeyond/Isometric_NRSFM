function [error_metric, T_poly, T_sel, T_norm] = test_two_view(dataset, pixel_noise, f, show_plot, decomp, choice, measure, use_visb, solver, grid, use_warp)

method = struct;

method.method = choice; 
% 0 for local approach: locally select the shape parameters
% 1 for global approach: select the solution by maximizing the consistency
% 2 for least median: select the least median (not a two-view method)

%%% local approach
method.measure = measure;
method.use_visb = use_visb; % flag for applying visibility condition
method.both_view = 0; % flag for penalizing on both views

%%% global approach
method.solver = solver;
method.c1 = 1.0; % parallelism for both views
method.c2 = 0.1; % least depth change in both views
method.c3 = 1.0; % visibility for both views

% parameters for constructing the Laplacian
method.sigma = 0.2;
method.ratio = 1.0; % threshold = ratio * average squared distance

dataset = ['./warps_data/', dataset];
if use_warp
    load(dataset);
else
    load(dataset, 'qgth', 'Pgth', 'Ngth', 'K');
end

%3D Ground truth points and normalized image points
n = length(qgth);
pairs_num = n - 1;
for i=1:n
    q_n(2*(i-1)+1:2*(i-1)+2,:) = qgth{i}(1:2, :);
end
%visiblity matrix: remove points visible in less than 3 views
visb = ones(n, size(q_n, 2));

%%% GROUND TRUTH NORMALS %%%%%
if exist('Ngth','var') == 0
    Ngth = create_gth_normals(Pgth,q_n,n, visb);
end

%%% PARAMETERS %%%%%
par = 2e-3; % schwarzian parameter.. needs to be tuned (usually its something close to 1e-3)
% %%%%% GRID OF POINTS %%%%
if grid %max(sum(visb)) == num
    % make a grid
    [I1u,I1v,I2u,I2v,visb] = create_grid(q_n,visb,20);
else
    % point-wise
    I1u = repmat(q_n(1,:),n-1,1);
    I1v = repmat(q_n(2,:),n-1,1);
    I2u = q_n(3:2:2*n,:);
    I2v = q_n(4:2:2*n,:);
end

% Add noise
if exist('K','var') == 0
    fx = f; fy = f;
else
    fx = K(1, 1); fy = K(2, 2);
end

% rng(2020);
num_p = length(I1u(1, :));
I1u(1, :) = randn(1, num_p) * pixel_noise / fx + I1u(1, :);
I1v(1, :) = randn(1, num_p) * pixel_noise / fy + I1v(1, :);
I2u = randn(pairs_num, num_p) * pixel_noise / fx + I2u;
I2v = randn(pairs_num, num_p) * pixel_noise / fy + I2v;


%%%%% SCHWARZIAN WARPS %%%%%
I1u = repmat(I1u(1, :), pairs_num, 1);
I1v = repmat(I1v(1, :), pairs_num, 1);
tic
if ~use_warp
    [I1u,I1v,I2u,I2v,J21a,J21b,J21c,J21d,J12a,J12b,J12c,J12d,H21uua,H21uub,H21uva,H21uvb,H21vva,H21vvb] = create_warps(I1u,I1v,I2u,I2v,visb,par, 1);
end
% [~,~,~,I1u,I1v,I2u,I2v,J21a,J21b,J21c,J21d,J12a,J12b,J12c,J12d,H21uua,H21uub,H21uva,H21uvb,H21vva,H21vvb] = create_tshirt_dataset(idx,length(idx), scene, par);
% create schwarzian warps for the dataset

toc

% view_id_list = 2:(1 + pairs_num);
view_id_list = [40];


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

method.eq_coef = eq_coef;
method.f1_coef = f1_coef;
method.f2_coef = f2_coef;
method.first_cubic_coef = first_cubic_coef;
method.t_all = t_all;
method.J21_all = J21_all;
method.visb = visb;
    
T_poly = zeros(1, length(view_id_list));
T_sel = zeros(1, length(view_id_list));
T_norm = zeros(1, length(view_id_list));

if decomp
    % coeff of k1       % coeff of k2        % constant term
    T1_12_k1 = -J21c;   T1_12_k2 = -J21d;    T1_12_c = (J12a.*H21uva + J12c.*H21uvb);%(H21vvb./J21d)/2;
    T2_12_k1 = -J21a;   T2_12_k2 = -J21b;    T2_12_c = (J12b.*H21uva + J12d.*H21uvb);%(H21uub./J21b)/2;

    % k1b = -T2_12 = a*k1 + b*k2 + t1;
    % k2b = -T1_12 = c*k1 + d*k2 + t2;
    a = -T2_12_k1; b = -T2_12_k2; c = -T1_12_k1; d = -T1_12_k2; t1 = -T2_12_c; t2 = -T1_12_c;

    [vec_W, vec_W_invt] = compute_W(repmat(I1u, pairs_num, 1), repmat(I1v, pairs_num, 1), I2u, I2v, a, b, c, d, t1, t2);
    [na, nb, na_, nb_] = compute_normal(vec_W, vec_W_invt);
end

count = 1;
% start solving all the cubics
for view_id = view_id_list
    idx = visb(1, :) == 1 & visb(view_id, :) == 1;
    
    method.view_id = view_id;
    I1u_tem = I1u(1, idx); I1v_tem = I1v(1, idx);
    I2u_tem = I2u(view_id - 1, idx); I2v_tem = I2v(view_id - 1, idx);
    t1 = t_all{view_id - 1}(1, idx);
    t2 = t_all{view_id - 1}(2, idx);

    Ngth_tem = [Ngth(1:3, :); Ngth((view_id * 3 - 2):(view_id * 3), :)];
    Pgth_tem = [Pgth(1:3, :); Pgth((view_id * 3 - 2):(view_id * 3), :)];
    q_n = [qgth{1}(1:2, :); qgth{view_id}(1:2, :)];

    tic
    if decomp
        num_p = length(t1);

        na_sub = reshape(na(view_id - 1, idx, :), num_p, 3)';
        nb_sub = reshape(nb(view_id - 1, idx, :), num_p, 3)';
        tem1 = I1u_tem .* na_sub(1, :) + I1v_tem .* na_sub(2, :) + na_sub(3, :);
        tem2 = I1u_tem .* nb_sub(1, :) + I1v_tem .* nb_sub(2, :) + nb_sub(3, :);
        k1 = [na_sub(1, :) ./ tem1;
            nb_sub(1, :) ./ tem2];
        k2 = [na_sub(2, :) ./ tem1;
            nb_sub(2, :) ./ tem2];

        na_sub_ = reshape(na_(view_id - 1, idx, :), num_p, 3)';
        nb_sub_ = reshape(nb_(view_id - 1, idx, :), num_p, 3)';
        tem1_ = I2u_tem .* na_sub_(1, :) + I2v_tem .* na_sub_(2, :) + na_sub_(3, :);
        tem2_ = I2u_tem .* nb_sub_(1, :) + I2v_tem .* nb_sub_(2, :) + nb_sub_(3, :);
        k1_ = [na_sub_(1, :) ./ tem1_;
            nb_sub_(1, :) ./ tem2_];
        k2_ = [na_sub_(2, :) ./ tem1_;
            nb_sub_(2, :) ./ tem2_];

        tem = nan(4, num_p);
        k1 = [tem; k1];
        k2 = [tem; k2];
        k1_ = [tem; k1_];
        k2_ = [tem; k2_];


    else
        % solve x1 and x2
        eq = eq_coef{view_id - 1}(:, idx);
        f1 = f1_coef{view_id - 1}(:, idx);
        f2 = f2_coef{view_id - 1}(:, idx);

        [x1, x2] = solve_cubic(eq, f1, f2, err);

        %%% Note that we replace all complex numbers with NaN
        %%% there are 6 roots, so x1 and x2 have 6 rows of solution

        J12a_all = repmat(J12a(view_id - 1, idx), 6, 1); J12b_all = repmat(J12b(view_id - 1, idx), 6, 1);
        J12c_all = repmat(J12c(view_id - 1, idx), 6, 1); J12d_all = repmat(J12d(view_id - 1, idx), 6, 1);

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
    mask = solution_selection(I1u_tem, I1v_tem, I2u_tem, I2v_tem, k1, k2, k1_, k2_, method);
    T_sel(count) = toc;

    % gather k1 and k2 for both views
    k1_all = [k1(mask)'; k1_(mask)'];
    k2_all = [k2(mask)'; k2_(mask)'];

    u_all = [I1u(1,idx);I2u_tem]; v_all = [I1v(1,idx);I2v_tem];

    % find normals on all surfaces N= [N1;N2;N3]
    N1 = k1_all; N2 = k2_all; N3 = 1-u_all.*k1_all-v_all.*k2_all;
    n = sqrt(N1.^2+N2.^2+N3.^2);
    N1 = N1./n ; N2 = N2./n; N3 = N3./n;

    N = [N1(:),N2(:),N3(:)]';
    N_res = reshape(N(:),3*num, size(u_all, 2));

    N_tem = zeros(3 * num, length(idx));
    N_tem(:, idx) = N_res;
    N_res = N_tem;
    
    tem_mask = ~isnan(k1);
    k1_all_another = [k1((~mask) & tem_mask)'; k1_((~mask) & tem_mask)'];
    k2_all_another = [k2((~mask) & tem_mask)'; k2_((~mask) & tem_mask)'];
    N1_ = k1_all_another; N2_ = k2_all_another; N3_ = 1-u_all.*k1_all_another-v_all.*k2_all_another;
    n = sqrt(N1_.^2+N2_.^2+N3_.^2);
    N1_ = N1_./n ; N2_ = N2_./n; N3_ = N3_./n;
    N_ = [N1_(:),N2_(:),N3_(:)]';
    N_res_another = reshape(N_(:),3*num, size(u_all, 2));
%     N_res = N_res_another;

    % compare GT normals
    er = 1e-5;
    t= 1e-3;
    nC = 40;
    
    umin = min(q_n(1,:))-t; umax = max(q_n(1,:))+t;
    vmin = min(q_n(2,:))-t; vmax = max(q_n(2,:))+t;
    
    bbs = bbs_create(umin, umax, nC, vmin, vmax, nC, 3);
    coloc = bbs_coloc(bbs, q_n(1,:), q_n(2,:));
    lambdas = er*ones(nC-3, nC-3);
    bending = bbs_bending(bbs, lambdas);
    cpts = (coloc'*coloc + bending) \ (coloc'*Pgth(1:3, :)');
    ctrlpts = cpts';
    
    dqu = bbs_eval(bbs, ctrlpts, I1u_tem',I1v_tem',1,0);
    dqv = bbs_eval(bbs, ctrlpts, I1u_tem',I1v_tem',0,1);
    
    nu = [dqu(1,:)./sqrt(sum(dqu.^2));dqu(2,:)./sqrt(sum(dqu.^2));dqu(3,:)./sqrt(sum(dqu.^2))];
    nv = [dqv(1,:)./sqrt(sum(dqv.^2));dqv(2,:)./sqrt(sum(dqv.^2));dqv(3,:)./sqrt(sum(dqv.^2))];
    nn = -cross(nu,nv);
    Ng1 = [nn(1,:)./sqrt(sum(nn.^2));nn(2,:)./sqrt(sum(nn.^2));nn(3,:)./sqrt(sum(nn.^2))];
      
    map1 = reshape(sqrt(sum(cross(Ng1(1:3, :), N_res(1:3, :)) .^ 2)), 20, 20);
    
    map1_another = reshape(sqrt(sum(cross(Ng1(1:3, :), N_res_another(1:3, :)) .^ 2)), 20, 20);
    figure
    histogram(map1(:))
    figure
    max(max(map1))
    map1 = map1 / max(max(map1));
    imshow(map1)
    
    
    umin = min(q_n(3,:))-t; umax = max(q_n(3,:))+t;
    vmin = min(q_n(4,:))-t; vmax = max(q_n(4,:))+t;
    
    bbs = bbs_create(umin, umax, nC, vmin, vmax, nC, 3);
    coloc = bbs_coloc(bbs, q_n(3,:), q_n(4,:));
    lambdas = er*ones(nC-3, nC-3);
    bending = bbs_bending(bbs, lambdas);
    cpts = (coloc'*coloc + bending) \ (coloc'*Pgth_tem(4:6, :)');
    ctrlpts = cpts';
    
    dqu = bbs_eval(bbs, ctrlpts, I2u_tem',I2v_tem',1,0);
    dqv = bbs_eval(bbs, ctrlpts, I2u_tem',I2v_tem',0,1);
    
    nu = [dqu(1,:)./sqrt(sum(dqu.^2));dqu(2,:)./sqrt(sum(dqu.^2));dqu(3,:)./sqrt(sum(dqu.^2))];
    nv = [dqv(1,:)./sqrt(sum(dqv.^2));dqv(2,:)./sqrt(sum(dqv.^2));dqv(3,:)./sqrt(sum(dqv.^2))];
    nn = -cross(nu,nv);
    Ng2 = [nn(1,:)./sqrt(sum(nn.^2));nn(2,:)./sqrt(sum(nn.^2));nn(3,:)./sqrt(sum(nn.^2))];
    map2 = reshape(sqrt(sum(cross(Ng2, N_res(4:6, :)) .^ 2)), 20, 20);
    map2_another = reshape(sqrt(sum(cross(Ng2, N_res_another(4:6, :)) .^ 2)), 20, 20);
    figure
    histogram(map2(:))    
    max(max(map2))
    map2 = map2 / max(max(map2));
    figure
    imshow(map2)
    
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




if length(view_id_list) > 1
    disp('Errors for the points')
    mean(1000 * error_metric(1:2, :)')
    disp('Errors for the normals (first and second)')
    mean(error_metric(3:4, :)')
    
    fprintf('Errors for the normals (all): %.3f \n', mean([mean(error_metric(3, :)), error_metric(4, :)]))
    fprintf('Errors for the depths (all): %.3f \n', 1000 * mean([mean(error_metric(1, :)), error_metric(2, :)]))
else
    disp('Errors for the points')
    disp(1000 * error_metric(1:2, :)')
    disp('Errors for the normals')
    error_metric(3:4, :)'
end
end
