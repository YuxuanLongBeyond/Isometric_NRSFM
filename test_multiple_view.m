function [err_n, err_p] = test_multiple_view(dataset, pixel_noise, f, solver, grid, show_plot, use_warp, views_num)

dataset = ['./warps_data/', dataset];
if use_warp
    load(dataset);
else
    load(dataset, 'qgth', 'Pgth', 'Ngth', 'K');
end

if ~strcmp(views_num, 'all')
    qgth = qgth(1:views_num);
    Pgth = Pgth(1:3*views_num, :);
end

pairs_num = length(qgth) - 1;
%3D Ground truth points and normalized image points
num = length(qgth);
for i=1:num
    q_n(2*(i-1)+1:2*(i-1)+2,:) = qgth{i}(1:2, :);
end
%visiblity matrix: remove points visible in less than 3 views
visb = ones(num, size(q_n, 2));

%%% GROUND TRUTH NORMALS %%%%%
if exist('Ngth','var') == 0
    Ngth = create_gth_normals(Pgth,q_n,num, visb);
end
if ~strcmp(views_num, 'all')
    Ngth = Ngth(1:3*views_num, :);
end

%%% PARAMETERS %%%%%
par = 2e-3; % schwarzian parameter.. needs to be tuned (usually its something close to 1e-3)
% 
% %%%%% GRID OF POINTS %%%%
if grid %max(sum(visb)) == num
    % make a grid
    [I1u,I1v,I2u,I2v,visb] = create_grid(q_n,visb,20);
else
    % point-wise
    I1u = repmat(q_n(1,:),num-1,1);
    I1v = repmat(q_n(2,:),num-1,1);
    I2u = q_n(3:2:2*num,:);
    I2v = q_n(4:2:2*num,:);
end
toc
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
I2u = randn(num - 1, num_p) * pixel_noise / fx + I2u;
I2v = randn(num - 1, num_p) * pixel_noise / fy + I2v;


%%%%% Compute warps %%%%%
I1u = repmat(I1u(1, :), pairs_num, 1);
I1v = repmat(I1v(1, :), pairs_num, 1);

if ~use_warp
    [I1u,I1v,I2u,I2v,J21a,J21b,J21c,J21d,J12a,J12b,J12c,J12d,H21uua,H21uub,H21uva,H21uvb,H21vva,H21vvb] = create_warps(I1u,I1v,I2u,I2v,visb,par, 0);
else
    if ~strcmp(views_num, 'all')
        J21a = J21a(1:(views_num - 1), :);
        J21b = J21b(1:(views_num - 1), :);
        J21c = J21c(1:(views_num - 1), :);
        J21d = J21d(1:(views_num - 1), :);
        J12a = J12a(1:(views_num - 1), :);
        J12b = J12b(1:(views_num - 1), :);
        J12c = J12c(1:(views_num - 1), :);
        J12d = J12d(1:(views_num - 1), :);
        H21uua = H21uua(1:(views_num - 1), :);
        H21uub = H21uub(1:(views_num - 1), :);
        H21uva = H21uva(1:(views_num - 1), :);
        H21uvb = H21uvb(1:(views_num - 1), :);
        H21vva = H21vva(1:(views_num - 1), :);
        H21vvb = H21vvb(1:(views_num - 1), :);
    end
end
disp('start computing normals...')
toc


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

if strcmp(solver, 'infP') || strcmp(solver, 'iso')
    e = 1+ I2u.^2 + I2v.^2; u = I2u; v = I2v;
    e1 = 1+ I1u.^2 + I1v.^2; u1 = I1u; v1 = I1v;

    % create sum of squares polynomial
    eq = create_polynomial_coefficients(a,b,c,d,t1,t2,e,e1,u,u1,v,v1);
    % minimise it to obtain depth derivatives
    res = solve_polynomial(eq);
    

%     % remove points with no solution
%     res_inf = res;
%     idx = find(res(:,1)==0); res(idx,:)=[];
%     u(:,idx) = []; v(:,idx) = []; u1(:,idx) = []; v1(:,idx) = []; e(:,idx)=[];e1(:,idx)=[];
%     a(:,idx) = []; b(:,idx) = []; c(:,idx) = []; d(:,idx) = []; t1(:,idx) = []; t2(:,idx) = [];
%     J12a(:,idx) = []; J12b(:,idx) = []; J12c(:,idx) = []; J12d(:,idx) = [];
%     H21uua(:,idx) = []; H21uub(:,idx) = []; H21uva(:,idx) = []; H21uvb(:,idx) = []; H21vva(:,idx) = []; H21vvb(:,idx) = [];

%     % recover first order derivatives on rest of the surfaces
%     k1_inf = [res(:,1)';a.*repmat(res(:,1)',num-1,1) + b.*repmat(res(:,2)',num-1,1) + t1];
%     k2_inf = [res(:,2)';c.*repmat(res(:,1)',num-1,1) + d.*repmat(res(:,2)',num-1,1) + t2];
%     k1_inf(:,idx)=[]; k2_inf(:,idx)=[];
%     k1_old = k1_inf; k2_old = k2_inf;
%     % NEW Christoffel symbols
%     % T1 = -2*k1 + k3*A;  T2 = k3*B;
%     % T3 = -k2 + k4*A;    T4 = -k1 + k4*B;
%     % T5 =  k5*A;         T6 = -2*k2 + k5*B;
%     res_old = res;

    if strcmp(solver, 'iso')
        disp('Solving general case: Alternate Approach.....')
        for iteration = 1:3
            iteration
            disp('Solve for second order derivatives.....')
            % second order derivatives k3 k4 k5
            [k3 k4 k5] = solve_k3k4k5(res_old,a,b,c,d,t1,t2,J12a,J12b,J12c,J12d,H21uua,H21uub,H21uva,H21uvb,H21vva,H21vvb,e,e1,u,u1,v,v1);
            k3(isnan(k3)) = 0; k4(isnan(k4)) = 0; k5(isnan(k5)) = 0;

            % recalculate first order derivatives using k3 k4 k5

            % find relation between first order derivatives of the surface i
            % with surface 1
            disp('recalculate first order derivatives.....')
            [a,b,c,d,t1,t2] = find_coefficient(k3,k4,k5,u,v,u1,v1,J21a,J21b,J21c,J21d,J12a,J12b,J12c,J12d,H21uua,H21uub,H21uva,H21uvb,H21vva,H21vvb);
            a(isnan(a)) = 0; b(isnan(b)) = 0; c(isnan(c)) = 0; d(isnan(d)) = 0; t1(isnan(t1)) = 0; t2(isnan(t2)) = 0;
            % k1_bar = a1*k1 + b1*k2 + t11;
            % k2_bar = c1*k1 + d1*k2 + t21;

            % solve for first order derivatives
            eq1 = create_polynomial_coefficients(a,b,c,d,t1,t2,e,e1,u,u1,v,v1);
            res = solve_polynomial(eq1);
            idx = find(res(:,1)==0); res(idx,:)= res_old(idx,:);
            res_new = res;
            k1_new = [res(:,1)';a.*repmat(res(:,1)',num-1,1) + b.*repmat(res(:,2)',num-1,1) + t1];
            k2_new = [res(:,2)';c.*repmat(res(:,1)',num-1,1) + d.*repmat(res(:,2)',num-1,1) + t2];

            th1=(k1_inf-k1_new).^2; th2=(k2_inf-k2_new).^2;
            blk = 1e-5;
            k1_final = k1_new.*(th1< blk) + k1_old.*(th1>= blk);
            k2_final = k2_new.*(th1< blk) + k2_old.*(th1>= blk);
            res_old = res;
            k1_old = k1_final; k2_old = k2_final;
        end

        k1_all = k1_final; k2_all = k2_final;

    else
        % recover first order derivatives on rest of the surfaces
        k1_all = [res(:,1)';a.*repmat(res(:,1)',num-1,1) + b.*repmat(res(:,2)',num-1,1) + t1];
        k2_all = [res(:,2)';c.*repmat(res(:,1)',num-1,1) + d.*repmat(res(:,2)',num-1,1) + t2];
    end
    
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
% 
%     % find indices with no solution
%     idx = find(res(:,1)==0);
%     N_res(:,idx) = []; u_all(:,idx) = []; v_all(:,idx) = [];

end

if strcmp(solver, 'fastDiffH')
    num_p = length(a(1, :));
    N_res = zeros(3 * num, num_p);
    [vec_W, vec_W_invt] = compute_W(I1u, I1v, I2u, I2v, a, b, c, d, t1, t2);
    [na, nb, ~, ~] = compute_normal(vec_W, vec_W_invt);
    
    err_tol = 0;
    
    for i = 1:num_p
        n1_all = reshape(na(:, i, :), pairs_num, 3)';
        n2_all = reshape(nb(:, i, :), pairs_num, 3)';
        x = [I1u(1, i);I1v(1, i);1];
        mask1 = (n1_all(3, :) ./ (x' * n1_all)) > -err_tol; 
        mask2 = (n2_all(3, :) ./ (x' * n2_all)) > -err_tol; 
        
        n_single = [n1_all(:, mask1 & (~mask2)), n2_all(:, (~mask1) & mask2)];
        P_single = ortho_projector(n_single);

        mask = mask1 & mask2;
        n1_all = n1_all(:, mask);
        n2_all = n2_all(:, mask);    
        N = cross(n1_all, n2_all)';
        [~, ~, V] = svd(N' * N); n = V(:, end);
     
%         w = 1;
%         for z = 1:10
%             N_ = sqrt(w) .* N;
%             [~, ~, V] = svd(N_' * N_); 
%             n = V(:, end);
%             w = 1 ./ abs(N * n);
%         end
        
%         N = [n1_all, n2_all];
%         [U, ~, ~] = svd(N * N'); n = U(:, 1); 
        
        iter_max = 1;
%         col_n = zeros(pairs_num, 3);
        for z = 1:iter_max
            dot1 = abs(n' * n1_all); dot2 = abs(n' * n2_all);
            mask = dot1 > dot2;
            candidate = [n1_all(:, mask), n2_all(:, ~mask)];
            P = [ortho_projector(candidate); P_single];

%             [U,S,V] = svd(P);
%             n = V(:, end);

            iter_max2 = 10;
            w = 1;
            for k = 1:iter_max2
                P_ = sqrt(w) .* P;
                [~,~,V] = svd(P_' * P_);
                n = V(:, end);
                tem = P * n;
                w = 1 ./ sqrt(tem(1:3:end) .^ 2 + tem(2:3:end) .^ 2 + tem(3:3:end) .^ 2);
                w = [w, w, w]';
                w = w(:);
                if any(isinf(w))
                    break
                end
            end
        end
%         n = mean(col_n)'; n = n / norm(n);
        N_res(1:3, i) = n;  

        n_(1, :) = vec_W{1, 1}(:, i) * n(1) + vec_W{1, 2}(:, i) * n(2) + vec_W{1, 3}(:, i) * n(3);
        n_(2, :) = vec_W{2, 1}(:, i) * n(1) + vec_W{2, 2}(:, i) * n(2) + vec_W{2, 3}(:, i) * n(3);
        n_(3, :) = vec_W{3, 1}(:, i) * n(1) + vec_W{3, 2}(:, i) * n(2) + vec_W{3, 3}(:, i) * n(3);
        
        n_ = n_ ./ sqrt(n_(1, :) .^ 2 + n_(2, :) .^ 2 + n_(3, :) .^ 2);
        N_res(4:end, i) = n_(:);
        
    end
    u_all = [I1u(1,:);I2u]; v_all = [I1v(1,:);I2v];
end


if strcmp(solver, 'polyH')
    num_p = length(a(1, :));
    N_res = zeros(3 * num, num_p);
    T = eye(3);
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
    u_all = [I1u(1,:);I2u]; v_all = [I1v(1,:);I2v];
end
toc


% Integrate normals to find depth
P_grid=calculate_depth(N_res,u_all,v_all,1e0);

% compare with ground truth
[P2,err_p] = compare_with_Pgth(P_grid,u_all,v_all,q_n,Pgth);
[N,err_n] = compare_with_Ngth(P2,q_n,Ngth);

% plot results
if show_plot
    for i=1:size(u_all,1)
         figure(i)
        plot3(Pgth(3*(i-1)+1,:),Pgth(3*(i-1)+2,:),Pgth(3*(i-1)+3,:),'go');
        hold on;
          plot3(P2(3*(i-1)+1,:),P2(3*(i-1)+2,:),P2(3*(i-1)+3,:),'ro');
          %quiver3(P2(3*(i-1)+1,:),P2(3*(i-1)+2,:),P2(3*(i-1)+3,:),N(3*(i-1)+1,:),N(3*(i-1)+2,:),N(3*(i-1)+3,:));
          %quiver3(Pgth(3*(i-1)+1,:),Pgth(3*(i-1)+2,:),Pgth(3*(i-1)+3,:),Ngth(3*(i-1)+1,:),Ngth(3*(i-1)+2,:),Ngth(3*(i-1)+3,:));
        hold off;
        axis equal;
    end
end

tem = mean(err_n');
fprintf('Average normal error on first frame: %.3f \n', mean(tem(1)))
fprintf('Average normal error on other frames: %.3f \n', mean(tem(2:end)))
fprintf('Average normal error on all frames: %.3f \n', mean(tem))


tem = mean(err_p');
fprintf('Average depth error on first frame: %.3f \n', mean(tem(1)) * 1000)
fprintf('Average depth error on other frames: %.3f \n', mean(tem(2:end)) * 1000)
fprintf('Average depth error on all frames: %.3f \n', mean(tem) * 1000)
end
