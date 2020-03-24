function mask = solution_selection(I1u, I1v, I2u, I2v, k1, k2, k1_, k2_, method)
n = length(I1u);
mask = zeros(6, n);
if method.method == 0
    measure = abs(k1).^2 + abs(k2).^2; % + abs(k1_).^2 + abs(k2_).^2;
    
    for i = 1:n
        tem = measure(:, i);
        tem = tem(~isnan(tem));
        mask(:, i) = abs(measure(:, i)) == min(abs(tem));
    end
      
end
if method.method == 1
    sigma = method.sigma;
    ratio = method.ratio;
    
%     F = [I1u; I1v; I2u; I2v];
% 
%     dist_map = zeros(n, n);
%     for i = 1:n
%         dist_map(:, i) = sum((F - F(:, i)) .^ 2);
%     end
% 
%     avg_dist = sum(sum(dist_map)) / (n ^ 2 - n);
% 
%     thre = ratio * avg_dist;
%     
%     dist_mask = dist_map > thre;
% 
%     W = exp(-dist_map / (2 * sigma ^ 2));
%     W = W - eye(n);
%     W(dist_mask) = 0;

    F1 = [I1u; I1v];
    F2 = [I2u; I2v];

    dist_map1 = zeros(n, n);
    dist_map2 = zeros(n, n);
    for i = 1:n
        dist_map1(:, i) = sum((F1 - F1(:, i)) .^ 2);
        dist_map2(:, i) = sum((F2 - F2(:, i)) .^ 2);
    end

    avg_dist1 = sum(sum(dist_map1)) / (n ^ 2 - n);
    thre1 = ratio * avg_dist1;
    dist_mask1 = dist_map1 > thre1;
    W1 = exp(-dist_map1 / (2 * sigma ^ 2));
    W1 = W1 - eye(n);
    W1(dist_mask1) = 0;
    
    avg_dist2 = sum(sum(dist_map2)) / (n ^ 2 - n);
    thre2 = ratio * avg_dist2;
    dist_mask2 = dist_map2 > thre2;
    W2 = exp(-dist_map2 / (2 * sigma ^ 2));
    W2 = W2 - eye(n);
    W2(dist_mask2) = 0;


%     % observe the non-sparsity
%     spy(W1)

%     fprintf('The average degree in the graph is %f \n', mean(sum(W1 ~= 0)))
%     fprintf('The average degree in the graph is %f \n', mean(sum(W2 ~= 0)))
    % graph Laplacian (unormalized form)
    L = diag(sum(W1)) - W1; 
    L_ = diag(sum(W2)) - W2;
    S1 = zeros(n, 6 * n); % allow maximum 6 real solutions
    S2 = S1; 
    S1_ = S1; S2_ = S1;

    index_list = zeros(1, n);
    for i = 1:n
        tem1 = k1(:, i);
        mask(:, i) = ~isnan(tem1);
        tem1 = tem1(~isnan(tem1));
        tem2 = k2(:, i);
        tem2 = tem2(~isnan(tem2));
        
        tem3 = k1_(:, i);
        tem3 = tem3(~isnan(tem3));
        
        tem4 = k2_(:, i);
        tem4 = tem4(~isnan(tem4));        

        index_list(i) = length(tem1);
        S1(i, (6 * i - 5):(6 * i - 6 + length(tem1))) = tem1;
        S2(i, (6 * i - 5):(6 * i - 6 + length(tem2))) = tem2;
        
        
        S1_(i, (6 * i - 5):(6 * i - 6 + length(tem3))) = tem3;
        S2_(i, (6 * i - 5):(6 * i - 6 + length(tem4))) = tem4;        

    end
    S1(:, sum(S1) == 0) = [];
    S2(:, sum(S2) == 0) = [];
    S1_(:, sum(S1_) == 0) = [];
    S2_(:, sum(S2_) == 0) = [];    

    B = S1; B(B ~= 0) = 1;
    C = diag(I1u) * S1 + diag(I1v) * S2;
    C_ = diag(I2u) * S1_ + diag(I2v) * S2_;
    A = S1' * L * S1 + S2' * L * S2 + C' * L * C; 
    % A = A + S1_' * L_ * S1_ + S2_' * L_ * S2_ + C_' * L_ * C_;


    m = size(A, 1);
    if strcmp(method.solver, 'qp1')
    %     invA = A ^ (-1);
    %     V = invA * B' * ((B * invA * B') \ ones(n, 1));
    
%     opts = struct;
%     opts.SYM = true;
%     tem = linsolve([A B'; B zeros(n, n)], [zeros(size(A, 1), 1);ones(n, 1)], opts);
    
        tem =  [A B'; B zeros(n, n)] \ [zeros(size(A, 1), 1);ones(n, 1)];
        V = tem(1:m);
    end
    
    if strcmp(method.solver, 'qp2')
        A = (A + A') / 2;
        V = quadprog(A, zeros(m, 1), -eye(m), zeros(m, 1), B, ones(n, 1));
    end


    if strcmp(method.solver, 'admm')
        tau = 0.001;
        lam = 1;
        max_iter = 10;
        mu = zeros(m, 1);
        W = zeros(m, 1);

        C = [A + tau * eye(m), B'; B, zeros(n, n)];
        C_inv = inv(C);
        C_inv_sub = C_inv(1:m, :);
        for i = 1:max_iter
            V = C_inv_sub * [tau * (W + mu); ones(n, 1)];

            z = V - mu;
            W = sign(z) .* max(abs(z) - lam / tau, 0);

            mu = mu + W - V;
        end        
    end
%     max_iter = 1000;

%     tau1 = 1; tau2 = 1;
%     alpha = 0.0001;
%     for i = 1:max_iter
%         tem_mask = V > 0.5;
%         V(tem_mask) = 1;
%         V(~tem_mask) = 0;        
%         grad = A * V;
%         V = V - alpha * grad;
%     end
%     k1 = S1 * V;
%     k2 = S2 * V;
%     k1_ = S1_ * V;
%     k2_ = S2_ * V;    

    cum_index = cumsum(index_list);
    j = 0;
    for i = 1:n
        subV = V((j + 1):(index_list(i) + j));
        j = cum_index(i);

        [~, id] = max(subV);
        tem = mask(:, i);
        find_index = find(tem == 1);
        tem(find_index(id)) = 2;
        mask(tem < 2, i) = 0;
    end
end

if method.method == 2
    pairs_num = length(method.eq_coef);
    
    
    for i = 1:n
        tem1 = k1(:, i);
        mask(:, i) = ~isnan(tem1);
        tem1 = tem1(~isnan(tem1));
        tem2 = k2(:, i);
        tem2 = tem2(~isnan(tem2));
        

        err_list = zeros(length(tem1), pairs_num);
        for j = 1:pairs_num
            if j ~= (method.view_id - 1)
                f1 = method.f1_coef{j};
                f2 = method.f2_coef{j};
                r = method.first_cubic_coef{j};
                J21_all = method.J21_all{j};
                
                J21 = J21_all(:, i); % a,b,c,d
                x1 = J21(1) * tem1 + J21(2) * tem2;
                x2 = J21(3) * tem1 + J21(4) * tem2;
                
                % first cubic (r0 + r1.*x1 + r2.*x2 + r3.*x1.*x2 + r4.*x1.^2 + r5.*x2.^2 + r6.*x1.*x2.^2 + r7.*x2.*x1.^2)
                % second cubic is f1 x1 + f2 = 0
                c1 = r(1, i) + r(2, i).*x1 + r(3, i).*x2 + r(4, i).*x1.*x2 + r(5, i).*x1.^2 + r(6, i).*x2.^2 + r(7, i).*x1.*x2.^2 + r(8, i).*x2.*x1.^2;
                c2 = polyval(flipud(f2(:, i))', x2) + x1 .* polyval(flipud(f1(:, i))', x2);
                err_list(:, j) = c1 .^ 2 + c2 .^ 2;
            end
        end
        err_list(:, method.view_id - 1) = [];
        
        [~, index] = min(median(err_list, 2));
        all_index = find(mask(:, i) == 1);
        mask(all_index(index), i) = 2;
        mask(:, i) = mask(:, i) == 2;
    end

end
mask = mask == 1;  
end

