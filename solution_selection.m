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
%     measure = max((k1 ./ k3) .^ 2 + (k2 ./ k3) .^ 2, (k1_ ./ k3_) .^ 2 + (k2_ ./ k3_) .^ 2);
    
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
    for i = 1:n
        tem1 = k1(:, i);
        mask(:, i) = ~isnan(tem1);
        tem1 = tem1(~isnan(tem1));
        tem2 = k2(:, i);
        tem2 = tem2(~isnan(tem2));
%         
%         tem3 = k1_(:, i);
%         tem3 = tem3(~isnan(tem3));
%         
%         tem4 = k2_(:, i);
%         tem4 = tem4(~isnan(tem4));        

        index_list(i) = length(tem1);
        
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
    r1 = k1 ./ k3;
    r2 = k2 ./ k3;
    
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
    
    S1 = cell(1, n);
    S2 = cell(1, n);
    for i = 1:n
        tem1 = r1(:, i);
        tem1 = tem1(~isnan(tem1));
        tem2 = r2(:, i);
        tem2 = tem2(~isnan(tem2));

        S1{i} = tem1;
        S2{i} = tem2;
    end    
    
    % pick the giant node as the seed
    [~, current_id] = max(sum(W1 ~= 0));
    
    labels = zeros(1, n); % label 0 if not selected
    all_index = 1:n;
    unique_solution_indicator = zeros(1, n);
    while any(labels == 0)
        r1_list = S1{current_id};
        r2_list = S2{current_id};
        
        index_visit = all_index(W1(current_id, :) > 0);
        
        num_choice = length(r1_list);
        if num_choice == 1
            unique_solution_indicator(current_id) = 1;
        end
        
        reward_matrix = zeros(num_choice, n);
        for i = 1:num_choice
            r1_tem = r1_list(i);
            r2_tem = r2_list(i);
            
            for j = index_visit
                r1_neighbour = S1{j};
                r2_neighbour = S2{j};
                reward = 0;
                if length(r1_neighbour) == 1
                    % more contribution for unique solution
                    reward = reward + W1(current_id, j) * exp(-abs(r1_tem - r1_neighbour)*100) * gain;
                    reward = reward + W1(current_id, j) * exp(-abs(r2_tem - r2_neighbour)*100) * gain;
                    unique_solution_indicator(j) = 1;
                else
                    reward = reward + W1(current_id, j) * exp(-min(abs(r1_tem - r1_neighbour))*100);
                    reward = reward + W1(current_id, j) * exp(-min(abs(r2_tem - r2_neighbour))*100);
                end
                reward_matrix(i, j) = reward;
            end
            
        end
        reward_total = sum(reward_matrix, 2);
        if sum(abs(mean(reward_total) - reward_total)) < 1e-3
            l1_norm = abs(r1_list) + abs(r2_list);
            [~, best_index] = min(l1_norm);
        else
            [~, best_index] = max(reward_total);
        end
        
        S1{current_id} = r1_list(best_index);
        S2{current_id} = r2_list(best_index);
        
        % zero the reward of unique solution
        reward_matrix(best_index, unique_solution_indicator == 1) = 0;
        if sum(reward_matrix(best_index, :)) == 0
            % special situation where surrounded by all unique solutions
            % so pick another one
            edge_num = sum(W1 ~= 0);
            edge_num(unique_solution_indicator == 1) = 0;
            if sum(edge_num) == 0
                break
            else
                [~, next_id] = max(edge_num);
            end
        else
            [~, next_id] = max(reward_matrix(best_index, :));
        end

        labels(current_id) = 1;
        current_id = next_id;
    end
    
    for i = 1:n
        r1_tem = S1{i};
        mask(:, i) = r1_tem == r1(:, i);
    end
    
end

mask = mask == 1;  
end

