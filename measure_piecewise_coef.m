function coef = measure_piecewise_coef(u, v, N)
n = length(u);

kNN = 10; % number of nearest neighbours
kNN_global = n - 1;

F = [u; v];

% compute distance maps for image coordinates
% and construct local neighbourhood graph and fully connected graph
W = zeros(n, n);
W_full = zeros(n, n);
dist_map = zeros(n, n);
for i = 1:n
    dist_map(:, i) = sum((F - F(:, i)) .^ 2);
    [~, index] = sort(dist_map(:, i));
    W(index(1:(kNN + 1)), i) = 1;
    W_full(index(1:(kNN_global + 1)), i) = 1;
    W(i, i) = 0;
    W_full(i, i) = 0;
    
    dot_diff(:, i) = acosd(N(:, i)' * N(:, W(:, i) == 1));
end

thre = median(median(dot_diff));

for i = 1:n
    mask = dot_diff(:, i) < thre;
    smoothness(i) = sum(2 - 2 * cosd(dot_diff(mask, i)));
    if sum(mask) > 0
        smoothness(i) = smoothness(i) / sum(mask);
    end
end
smoothness(smoothness == 0) = max(smoothness);

L_full = diag(sum(W_full)) - W_full;
global_coef = trace(N * L_full * N') / n / kNN_global;

% coef = exp(100 * local_coef / global_coef);
% coef = 100 * (smoothness / global_coef);
coef = (500 * smoothness / global_coef);
    
end