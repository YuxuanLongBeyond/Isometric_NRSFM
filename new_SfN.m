function [ctrpts3Dn]=new_SfN(bbs,coloc,bending,m1,n1, weight)
coloc_du = bbs_coloc_deriv(bbs, m1(1,:)', m1(2,:)', 1, 0);
coloc_dv = bbs_coloc_deriv(bbs, m1(1,:)', m1(2,:)', 0, 1);
dot = sum(n1 .* m1)';

A = [coloc_du .* dot + n1(1, :)' .* coloc; coloc_dv .* dot + n1(2, :)' .* coloc];

% design the weight matrix
% n1 = n1';
% W = dot ./ sqrt(n1(:, 1) .^ 2 + n1(:, 2) .^ 2 + 1e-6);
% W = (n1(:, 3) .^ 2) ./ (n1(:, 1) .^ 2 + n1(:, 2) .^ 2 + 1e-3);
W = weight';
% W = dot;

A = A .* [W; W];

% % 
% [m, n] = size(A);
% f = [zeros(n, 1); ones(m, 1)];
% H = zeros(m + n, m + n); H(1:n, 1:n) = 2 * bending;
% A_ineq = [A, -eye(m); -A, -eye(m)]; b_ineq = zeros(2 * m, 1);
% A_eq = [ones(1, n), zeros(1, m)]; b_eq = 1;
% c = quadprog(H, f, A_ineq, b_ineq, A_eq, b_eq);
% ctrpts3Dn =  c(1:n)';

[~, n] = size(A);
H = bending + A' * A;
A_eq = ones(1, n); b_eq = 1;
c = quadprog(H, [], [], [], A_eq, b_eq);
ctrpts3Dn =  c';

end