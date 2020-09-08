function P = ortho_projector(n_all)
num = size(n_all, 2);
P = zeros(3 * num, 3);
for k = 1:num
    n = n_all(:, k);
    P((3 * (k - 1) + 1):(3 * k), :) = eye(3) - n * n';
end
end