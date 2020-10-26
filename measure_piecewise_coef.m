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
    
    
    smoothness(i) = norm((N(:, W(:, i) == 1) - N(:, i)), 'fro') ^ 2 / kNN;
    
end

L_full = diag(sum(W_full)) - W_full;
global_coef = trace(N * L_full * N') / n / kNN_global;

% coef = exp(100 * local_coef / global_coef);
coef = 100 * (smoothness / global_coef);
    
end