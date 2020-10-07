function coef = measure_smoothness(u, v, N)
n = length(u);

kNN = 10;
kNN_global = n - 1;

F = [u; v];

% compute distance maps for image coordinates
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
end

L = diag(sum(W)) - W;
local_coef = trace(N * L * N') / n / kNN;

L_full = diag(sum(W_full)) - W_full;
global_coef = trace(N * L_full * N') / n / kNN_global;

% W_full = ones(n, n) - eye(n);
% global_coef = trace(N * (diag(sum(W_full)) - W_full) * N') / n / (n - 1);

% coef = exp(100 * local_coef / global_coef);

coef = exp(10 * local_coef / global_coef);
    
end