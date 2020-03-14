function mask = solution_selection(I1u, I1v, I2u, I2v, k1, k2, method)
n = length(I1u);
mask = zeros(6, n);
if method.method == 0
    measure = abs(k1).^2 + abs(k2).^2;
    
    for i = 1:n
        tem = measure(:, i);
        tem = tem(~isnan(tem));
        mask(:, i) = abs(measure(:, i)) == min(abs(tem));
    end
      
end

if method.method == 1
    sigma = method.sigma;
    ratio = method.ratio;
    tao = method.tao;
    
    F = [I1u; I1v; I2u; I2v];

    dist_map = zeros(n, n);
    for i = 1:n
        dist_map(:, i) = sum((F - F(:, i)) .^ 2);
    end

    avg_dist = sum(sum(dist_map)) / (n ^ 2 - n);
    thre = ratio * avg_dist;
    
    dist_mask = dist_map > thre;

    W = exp(-dist_map / (2 * sigma ^ 2));
    W = W - eye(n);
    W(dist_mask) = 0;

    
    
    % observe the non-sparsity
    spy(W)

    fprintf('The average degree in the graph is %f \n', mean(sum(W ~= 0)))

    % graph Laplacian (unormalized form)
    L = diag(sum(W)) - W;
    S1 = zeros(n, 3 * n); % allow maximum 3 real solutions
    S2 = S1;
    index_list = zeros(1, n);
    for i = 1:n
        tem1 = k1(:, i);
        mask(:, i) = ~isnan(tem1);
        tem1 = tem1(~isnan(tem1));
        tem2 = k2(:, i);
        tem2 = tem2(~isnan(tem2));

        index_list(i) = length(tem1);
        S1(i, (3 * i - 2):(3 * i - 3 + length(tem1))) = tem1;
        S2(i, (3 * i - 2):(3 * i - 3 + length(tem2))) = tem2;

    end

    S1(:, sum(S1) == 0) = [];
    S2(:, sum(S2) == 0) = [];

    B = S1; B(B > 0) = 1;
    C = diag(I1u) * S1 + diag(I1v) * S2; 
    A = S1' * L * S1 + S2' * L * S2 + C' * L * C;
    V =  ((tao * (B)' * B + A) \ (B' * ones(n, 1))) * tao;

%     lambda = 0.1;
%     V = ones(size(A, 1), 1) * Inf;
%     for i = 1:10
%         tem = [A B'; B zeros(n, n)] \ [lambda ./ V;ones(n, 1)];
%         V = tem(1:size(A, 1));
%     end
    
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

mask = mask == 1;  
end
