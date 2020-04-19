function mask = solution_selection(I1u, I1v, I2u, I2v, k1, k2, k1_, k2_, method)
n = length(I1u);
mask = zeros(6, n); % mask (boolean array) should indicate the most desirable real solution
if method.method == 0
    % l2 norm
%     measure = abs(k1).^2 + abs(k2).^2; % + abs(k1_).^2 + abs(k2_).^2;
% 
    I1u_all = repmat(I1u, 6, 1); I1v_all = repmat(I1v, 6, 1);
    I2u_all = repmat(I2u, 6, 1); I2v_all = repmat(I2v, 6, 1);

    
    k3 = 1 - I1u_all .* k1 - I1v_all .* k2;
    k3_ = 1 - I2u_all .* k1_ - I2v_all .* k2_;
    measure = (k1 ./ k3) .^ 2 + (k2 ./ k3) .^ 2;
%     measure = measure + (k1_ ./ k3_) .^ 2 + (k2_ ./ k3_) .^ 2;
    
    lam = 0.0;
    measure = measure + lam * (abs(k1_ ./ k3_ - k1 ./ k3) + abs(k2_ ./ k3_ - k2 ./ k3));
%     measure = k1 .^ 2 + k2 .^ 2 + (1 - I1u_all .* k1 - I1v_all .* k2) .^ 2;
    for i = 1:n
        tem = measure(:, i);
        tem = tem(~isnan(tem)); % only select the real roots
        % (note we have set the complex numbers as NaN)
        % we choose the one with least norm
        mask(:, i) = abs(measure(:, i)) == min(abs(tem));
    end
    
end
if method.method == 1
    sigma = method.sigma;
    ratio = method.ratio;

    F1 = [I1u; I1v]; % image coordinates on first view
%     F2 = [I2u; I2v]; % image coordinates on second view

    % compute distance maps for image coordinates
    dist_map1 = zeros(n, n);
%     dist_map2 = zeros(n, n);
    for i = 1:n
        dist_map1(:, i) = sum((F1 - F1(:, i)) .^ 2);
%         dist_map2(:, i) = sum((F2 - F2(:, i)) .^ 2);
    end

    % compute the adjacency matrix for first view
    avg_dist1 = sum(sum(dist_map1)) / (n ^ 2 - n);
    thre1 = ratio * avg_dist1;
    dist_mask1 = dist_map1 > thre1;
    W1 = exp(-dist_map1 / (2 * sigma ^ 2));
    W1 = W1 - eye(n); % leave zero on diagonal
    W1(dist_mask1) = 0;
    
%     % compute the adjacency matrix for second view
%     avg_dist2 = sum(sum(dist_map2)) / (n ^ 2 - n);
%     thre2 = ratio * avg_dist2;
%     dist_mask2 = dist_map2 > thre2;
%     W2 = exp(-dist_map2 / (2 * sigma ^ 2));
%     W2 = W2 - eye(n);
%     W2(dist_mask2) = 0;

%     spy(W1)

    % compute the graph Laplacians (unormalized form)
    L = diag(sum(W1)) - W1; 
%     L_ = diag(sum(W2)) - W2;

    S1 = zeros(n, 6 * n); % allow maximum 6 real solutions
    S2 = S1; 
%     S1_ = S1; S2_ = S1;
    start_index = 1;
    index_list = zeros(1, n);
    template_index = 1:6;
    non_nan_index = cell(1, n);
    for i = 1:n
        tem1 = k1(:, i);
        tem_mask = ~isnan(tem1);
        mask(:, i) = tem_mask;
        tem1 = tem1(tem_mask);
        tem2 = k2(:, i);
        tem2 = tem2(tem_mask);
%         
%         tem3 = k1_(:, i);
%         tem3 = tem3(~isnan(tem3));
%         
%         tem4 = k2_(:, i);
%         tem4 = tem4(~isnan(tem4));        

        index_list(i) = length(tem1);
        non_nan_index{i} = template_index(tem_mask);
        % collect the real solutions into the row
        S1(i, start_index:(start_index + length(tem1) - 1)) = tem1;
        S2(i, start_index:(start_index + length(tem2) - 1)) = tem2; 
        start_index = start_index + index_list(i);

    end
    S1 = S1(:, 1:(start_index - 1));
    S2 = S2(:, 1:(start_index - 1));

    B = S1; B(B ~= 0) = 1;
    C = diag(I1u) * S1 + diag(I1v) * S2;
%     C_ = diag(I2u) * S1_ + diag(I2v) * S2_;
    A = S1' * L * S1 + S2' * L * S2 + C' * L * C;
    
    %%% currently we only use the graph information from first view
    % A = A + S1_' * L_ * S1_ + S2_' * L_ * S2_ + C_' * L_ * C_;

    m = size(A, 1);
    if strcmp(method.solver, 'irqp')
        opts = struct;
        opts.SYM = true;        
        
        max_iter = 10;
        V = zeros(m, 1);
        for i = 1:max_iter
            invd = 1 ./ (ones(n, 1) - C * V);
            L_ = invd .* L .* (invd');
            A = S1' * L_ * S1 + S2' * L_ * S2;
%             tem = linsolve([A B'; B zeros(n, n)], [zeros(size(A, 1), 1);ones(n, 1)], opts);
%             V = tem(1:m);
            
            A = (A + A') / 2;
            V = quadprog(A, zeros(m, 1), -eye(m), zeros(m, 1), B, ones(n, 1));            
        end
    end

    if strcmp(method.solver, 'qp1')
        % no inequality constraint

        opts = struct;
        opts.SYM = true;
        
        % solve KKT equation
        tem = linsolve([A B'; B zeros(n, n)], [zeros(size(A, 1), 1);ones(n, 1)], opts);

%         tem =  [A B'; B zeros(n, n)] \ [zeros(size(A, 1), 1);ones(n, 1)];
        V = tem(1:m); % extract V since tem also consists dual variable
    end

    
    if strcmp(method.solver, 'qp2')
        % solve standard QP
        A = (A + A') / 2;
        V = quadprog(A, zeros(m, 1), -eye(m), zeros(m, 1), B, ones(n, 1));
    end


    if strcmp(method.solver, 'admm')
        % use L1 norm to replace the inequality constraint
        
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
    % once V is obtained, we fill out the mask
    j = 1;
    for i = 1:n
        subV = V(j:(j + index_list(i) - 1)); % v_i
        j = j + index_list(i);
        [~, id] = max(subV); % the maximum of v_i indicates the most desirable solution
        tem = mask(:, i);
        find_index = find(tem == 1);
        tem(find_index(id)) = 2;
        mask(tem < 2, i) = 0;
    end
end

if method.method == 2
    %%% least median method
    
    pairs_num = length(method.eq_coef);
    
    
    for i = 1:n
        
        % extract the real roots
        tem1 = k1(:, i);
        mask(:, i) = ~isnan(tem1);
        tem1 = tem1(~isnan(tem1));
        tem2 = k2(:, i);
        tem2 = tem2(~isnan(tem2));
       
        % ready to collect the sum of squared errors
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
                
                % first cubic
                c1 = r(1, i) + r(2, i).*x1 + r(3, i).*x2 + r(4, i).*x1.*x2 + r(5, i).*x1.^2 + r(6, i).*x2.^2 + r(7, i).*x1.*x2.^2 + r(8, i).*x2.*x1.^2;
                
                % second cubic
                c2 = polyval(flipud(f2(:, i))', x2) + x1 .* polyval(flipud(f1(:, i))', x2);
                err_list(:, j) = c1 .^ 2 + c2 .^ 2; % sum of squared errors
            end
        end
        err_list(:, method.view_id - 1) = [];
        
        [~, index] = min(median(err_list, 2)); % choose least median
        
        all_index = find(mask(:, i) == 1);
        mask(all_index(index), i) = 2;
        mask(:, i) = mask(:, i) == 2;
    end

end


if method.method == 3

    I1u_all = repmat(I1u, 6, 1); I1v_all = repmat(I1v, 6, 1);
    k3 = 1 - I1u_all .* k1 - I1v_all .* k2;    
    
    gain = 1.05;
    
    sigma = 0.1;
    ratio = 0.1;
    

    F1 = [I1u; I1v]; % image coordinates on first view
%     F2 = [I2u; I2v]; % image coordinates on second view

    % compute distance maps for image coordinates
    dist_map1 = zeros(n, n);
%     dist_map2 = zeros(n, n);
    for i = 1:n
        dist_map1(:, i) = sum((F1 - F1(:, i)) .^ 2);
%         dist_map2(:, i) = sum((F2 - F2(:, i)) .^ 2);
    end

    % compute the adjacency matrix for first view
    avg_dist1 = sum(sum(dist_map1)) / (n ^ 2 - n);
    thre1 = ratio * avg_dist1;
    dist_mask1 = dist_map1 > thre1;
    W1 = exp(-dist_map1 / (2 * sigma ^ 2));
    W1 = W1 - eye(n); % leave zero on diagonal
    W1(dist_mask1) = 0;   
%     spy(W1)
    
    template_index = 1:6;
    non_nan_index = cell(1, n);
    S = cell(1, n);
    for i = 1:n
        tem1 = k1(:, i);
        tem_mask = ~isnan(tem1);
        tem1 = tem1(tem_mask);
        tem2 = k2(:, i);
        tem2 = tem2(tem_mask);
        tem3 = k3(:, i);
        tem3 = tem3(tem_mask);        

        non_nan_index{i} = template_index(tem_mask);
        
        normal = [tem1, tem2, tem3];
        S{i} = normal ./ sqrt(sum(normal .^ 2, 2));
%         S{i} = normal;
    end    
%     
%     % pick the giant node as the seed
%     [~, current_id] = max(sum(W1 ~= 0));
%     
%     labels = zeros(1, n); % label 0 if not selected
%     unique_solution_indicator = zeros(1, n);
% %     count = 0;
%     while any(labels == 0)
% %         if count == 10
% %             break
% %         end
% %         count = count + 1;
%         
%         n1 = S{current_id};
%         num_choice = size(n1, 1);
%         if num_choice == 1
%             unique_solution_indicator(current_id) = 1;
%         end
%         
%         index_visit = find(W1(current_id, :) > 0);
%         reward_matrix = zeros(num_choice, n);
%         for j = index_visit
%             n2 = S{j};
%             reward = 0;
% %             r1_diff = abs(n1(:, 1) ./ n1(:, 3) - (n2(:, 1) ./ n2(:, 3))');
% %             r2_diff = abs(n1(:, 2) ./ n1(:, 3) - (n2(:, 2) ./ n2(:, 3))');
%             if size(n2, 1) == 1
% %                 reward = reward + W1(current_id, j) * exp(-(r1_diff + r2_diff) * 100) * gain;
% %                 reward = reward + W1(current_id, j) * exp(-r2_diff * 100) * gain;
%                 reward = reward + W1(current_id, j) * exp((2 * n1 * n2' - 2) * 10) * gain;
%                 unique_solution_indicator(j) = 1;
%             else
% %                 reward = reward + W1(current_id, j) * exp(-min(r1_diff + r2_diff, [], 2) * 100);
% %                 reward = reward + W1(current_id, j) * exp(-min(r2_diff, [], 2) * 100);
%                 reward = reward + W1(current_id, j) * exp((2 * max(n1 * n2', [], 2) - 2) * 10);
%             end
%             reward_matrix(:, j) = reward;
%         end
%         
% 
%         
%         if labels(current_id) == 0
%             reward_total = sum(reward_matrix, 2);
%             if sum(abs(mean(reward_total) - reward_total)) < 1e-3
%                 measure = (n1(:, 1) ./ n1(:, 3)) .^ 2 + (n1(:, 2) ./ n1(:, 3)) .^ 2;
%                 [~, best_index] = min(measure);
%             else
%                 [~, best_index] = max(reward_total);
%             end        
%             % update            
%             
%             S{current_id} = n1(best_index, :);
%             mask(non_nan_index{current_id}(best_index), current_id) = 2;
%             mask(:, current_id) = mask(:, current_id) > 1;    
%         else
%             best_index = 1;
%         end
% 
%         
%         % zero the reward of unique solution
%         reward_matrix(best_index, unique_solution_indicator == 1) = 0;
%         if sum(reward_matrix(best_index, :)) == 0
%             % special situation where surrounded by all unique solutions
%             % so pick another one
%             edge_num = sum(W1 ~= 0);
%             edge_num(unique_solution_indicator == 1) = 0;
%             if sum(edge_num) == 0
%                 break
%             else
%                 [~, next_id] = max(edge_num);
%             end
%         else
%             [~, next_id] = max(reward_matrix(best_index, :));
%         end
% 
%         labels(current_id) = 1;
%         current_id = next_id;
%         
%     end
% 
%     
    
    [~, ind_list] = sort(sum(W1 > 0), 'descend');
    for ind = ind_list
        n1 = S{ind};
        if size(n1, 1) > 1
            index_visit = find(W1(ind, :) > 0);
            reward = zeros(size(n1, 1), 1);
            for j = index_visit
                n2 = S{j};
                if size(n2, 1) == 1
                    reward = reward + W1(ind, j) * ((2 * n1 * n2' - 2) * 1) * gain;
                else
                    reward = reward + W1(ind, j) * ((2 * max(n1 * n2', [], 2) - 2) * 1);
                end
            end
            if sum(abs((mean(reward) - reward))) < 1e-3
                measure = (n1(:, 1) ./ n1(:, 3)) .^ 2 + (n1(:, 2) ./ n1(:, 3)) .^ 2;
                [~, max_index] = min(measure);
            else
                [~, max_index] = max(reward);
            end
            S{ind} = n1(max_index, :);
        end
        mask(non_nan_index{ind}(max_index), ind) = 2;
        mask(:, ind) = mask(:, ind) > 1;
    end
end

if method.method == 4
    I1u_all = repmat(I1u, 6, 1); I1v_all = repmat(I1v, 6, 1);
    I2u_all = repmat(I2u, 6, 1); I2v_all = repmat(I2v, 6, 1);
    k3 = 1 - I1u_all .* k1 - I1v_all .* k2;    
    k3_ = 1 - I2u_all .* k1_ - I2v_all .* k2_;
    sigma = method.sigma;
    ratio = method.ratio;
    

    F1 = [I1u; I1v]; % image coordinates on first view
    F2 = [I2u; I2v]; % image coordinates on second view

    % compute distance maps for image coordinates
    dist_map1 = zeros(n, n);
    dist_map2 = zeros(n, n);
    for i = 1:n
        dist_map1(:, i) = sum((F1 - F1(:, i)) .^ 2);
        dist_map2(:, i) = sum((F2 - F2(:, i)) .^ 2);
    end

    % compute the adjacency matrix for first view
    avg_dist1 = sum(sum(dist_map1)) / (n ^ 2 - n);
    thre1 = ratio * avg_dist1;
    dist_mask1 = dist_map1 > thre1;
    W1 = exp(-dist_map1 / (2 * sigma ^ 2));
    W1 = W1 - eye(n); % leave zero on diagonal
    W1(dist_mask1) = 0;
    
    % compute the adjacency matrix for second view
    avg_dist2 = sum(sum(dist_map2)) / (n ^ 2 - n);
    thre2 = ratio * avg_dist2;
    dist_mask2 = dist_map2 > thre2;
    W2 = exp(-dist_map2 / (2 * sigma ^ 2));
    W2 = W2 - eye(n);
    W2(dist_mask2) = 0;

%     spy(W1)

    % compute the graph Laplacians (unormalized form)
    L = diag(sum(W1)) - W1; 
    L_ = diag(sum(W2)) - W2;

    
    template_index = 1:6;
    non_nan_index = cell(1, n);
    S1 = zeros(n, 6 * n); % allow maximum 6 real solutions
    S2 = S1; S3 = S1;
    S1_ = zeros(n, 6 * n); % allow maximum 6 real solutions
    S2_ = S1_; S3_ = S1_;    
    R1 = S1; R2 = R1;
    R1_ = R1; R2_ = R1;
    
    start_index = 1;
    index_list = zeros(1, n);    
    measure = (k1 ./ k3) .^ 2 + (k2 ./ k3) .^ 2;
    V0 = zeros(6 * n, 1);
    for i = 1:n
        tem1 = k1(:, i);
        tem_mask = ~isnan(tem1);
        tem1 = tem1(tem_mask);
        tem2 = k2(:, i);
        tem2 = tem2(tem_mask);
        tem3 = k3(:, i);
        tem3 = tem3(tem_mask);        
        
        tem1_ = k1_(:, i);
        tem1_ = tem1_(tem_mask);
        tem2_ = k2_(:, i);
        tem2_ = tem2_(tem_mask);
        tem3_ = k3_(:, i);
        tem3_ = tem3_(tem_mask);              

        num_sol = length(tem1);
        [~, id] = min(measure(tem_mask, i));
        V0(start_index + id - 1) = 1;
        
        index_list(i) = num_sol;
        
        non_nan_index{i} = template_index(tem_mask);
        
        normal = [tem1, tem2, tem3];
        normal = normal ./ sqrt(sum(normal .^ 2, 2));
        
        normal_ = [tem1_, tem2_, tem3_];
        normal_ = normal_ ./ sqrt(sum(normal_ .^ 2, 2));

        
        S1(i, start_index:(start_index + num_sol - 1)) = normal(:, 1);
        S2(i, start_index:(start_index + num_sol - 1)) = normal(:, 2);
        S3(i, start_index:(start_index + num_sol - 1)) = normal(:, 3);
        R1(i, start_index:(start_index + num_sol - 1)) = normal(:, 1) ./ normal(:, 3);
        R2(i, start_index:(start_index + num_sol - 1)) = normal(:, 2) ./ normal(:, 3);
        
        
        S1_(i, start_index:(start_index + num_sol - 1)) = normal_(:, 1);
        S2_(i, start_index:(start_index + num_sol - 1)) = normal_(:, 2);
        S3_(i, start_index:(start_index + num_sol - 1)) = normal_(:, 3);
        R1_(i, start_index:(start_index + num_sol - 1)) = normal_(:, 1) ./ normal_(:, 3);
        R2_(i, start_index:(start_index + num_sol - 1)) = normal_(:, 2) ./ normal_(:, 3);        
        
        start_index = start_index + index_list(i);        
    end    

    S1 = S1(:, 1:(start_index - 1));
    S2 = S2(:, 1:(start_index - 1));    
    S3 = S3(:, 1:(start_index - 1));    
    S1_ = S1_(:, 1:(start_index - 1));
    S2_ = S2_(:, 1:(start_index - 1));    
    S3_ = S3_(:, 1:(start_index - 1));        
    
    R1 = R1(:, 1:(start_index - 1));    
    R2 = R2(:, 1:(start_index - 1));      
    R1_ = R1_(:, 1:(start_index - 1));    
    R2_ = R2_(:, 1:(start_index - 1));      
    
    
    V0 = V0(1:(start_index - 1));   
    
    A_original = S1' * L * S1 + S2' * L * S2 + S3' * L * S3;
    loss = max(V0' * A_original * V0, 0);
    c0 = 1.0;
    c1 = 1.0;
    c2 = 1.0;
    c3 = 1.0;
    
    B = S1; B(B ~= 0) = 1;
    A = c0 * A_original;
    A = A + c2 * (S1_' * L_ * S1_ + S2_' * L_ * S2_ + S3_' * L_ * S3_);
    A = A + c1 * (R1' * R1 + R2' * R2); % regularization
    A = A + c3 * (R1_' * R1_ + R2_' * R2_);
    m = size(A, 1);
    A = A / n;
%     spy(A)
    % use L1 norm to replace the inequality constraint
    s = 0.1;
    tau = 5;
    lam1 = 10;
    lam2 = lam1 * s;
    max_iter = 40;
    mu = zeros(m, 1);
    W = V0;

    A = A - eye(m) * lam2;

    C = [A + tau * eye(m), B'; B, zeros(n, n)];
    C_inv = inv(C);
    C_inv_sub = C_inv(1:m, :);
    if loss > 0
        for i = 1:max_iter
            V = C_inv_sub * [tau * (W + mu); ones(n, 1)];


%             if i == 1
%                 loss_new = V' * A_original * V;
%                 if loss_new >= loss
%                     V = V0;
%                     break
%                 end
%             end

            z = V - mu;
            W = sign(z) .* max(abs(z) - lam1 / tau, 0);

            mu = mu + W - V;
        end
        
        j = 1;
        for i = 1:n
            subV = V(j:(j + index_list(i) - 1)); % v_i
            tem_mask = max(subV) == subV;
            subV(tem_mask) = 1;
            subV(~tem_mask) = 0;
            V(j:(j + index_list(i) - 1)) = subV;
            j = j + index_list(i);
        end
        loss_new = V' * A_original * V;
        if loss_new >= loss
            V = V0;
        end

    else
        V = V0;
    end

    j = 1;
    for i = 1:n
        subV = V(j:(j + index_list(i) - 1)); % v_i
%         if all(abs(mean(subV) - subV) < 1e-3)
%             subV = V0(j:(j + index_list(i) - 1));
%         end
        j = j + index_list(i);
        mask(:, i) = zeros(6, 1);
        mask(non_nan_index{i}(subV == 1), i) = 1;
        
    end
    
end


mask = mask == 1;  
end

