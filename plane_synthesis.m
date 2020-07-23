clear all
close all

% make up an ideal camera
f = 1000;
cx = 0;
cy = 0;

sample_num = 2; % number of samples per row or column in the image
range_start = -200;
range_end = 200;
pairs_num = 2;

% define normal such that normal^T X = a1 where a1 = 1
n1 = [0, 0, 1]';
% n1 = [0.8, 0, 0.2]';
% n1 = [0.3, 0.3, 0.8]';
% n1 = [0.2 0.2 0.8]';
% n1 = [0.25 0.25 0.8]';
% n1 = [0.4 -0.6 0.8]';
n1 = n1 / norm(n1);
d1 = 10;
n1 = n1 / d1;

total_num = sample_num ^ 2;

% create image points in the first view
[x1, y1] = meshgrid(linspace(range_start, range_end, sample_num), linspace(range_start, range_end, sample_num));
x1 = reshape(x1, 1, total_num);
y1 = reshape(y1, 1, total_num);

K = [f, 0, cx; 0, f, cy; 0, 0, 1];

q2 = K \ [x1; y1; ones(1, total_num)];

R = rodrigues(randn(1, 3) * 0.2);
n1 = R * n1;
disp(n1 / norm(n1))

% q2 = R * (d1 * q2);
q2 = [0, n1(3), -n1(2), 0;-n1(3),0, n1(1), 0;n1(2),-n1(1), 1, n1(2)];

q2(1, :) = q2(1, :) ./ q2(3, :);
q2(2, :) = q2(2, :) ./ q2(3, :);
I1u = q2(1, :);
I1v = q2(2, :);



% create a set of random rotation, translation
t_all = randn(pairs_num, 3) * 0.5;
% t_all = randn(pairs_num, 3) * 0.0;
r_all = randn(pairs_num, 3) * 0.5;

% create image points in the second view
q1 = [I1u; I1v; ones(1, total_num)];
qgth = cell(1, pairs_num + 1);
qgth{1} = q1;
R_set = cell(1, pairs_num);
I2u = zeros(pairs_num, total_num);
I2v = zeros(pairs_num, total_num);
H21uua = zeros(pairs_num, total_num);
H21uub = zeros(pairs_num, total_num);
H21uva = zeros(pairs_num, total_num);
H21uvb = zeros(pairs_num, total_num);
H21vva = zeros(pairs_num, total_num);
H21vvb = zeros(pairs_num, total_num);
J12a = zeros(pairs_num, total_num);
J12b = zeros(pairs_num, total_num);
J12c = zeros(pairs_num, total_num);
J12d = zeros(pairs_num, total_num);
J21a = zeros(pairs_num, total_num);
J21b = zeros(pairs_num, total_num);
J21c = zeros(pairs_num, total_num);
J21d = zeros(pairs_num, total_num);
Ngth = zeros(3 * (pairs_num + 1), total_num);
Ngth(1:3, :) = -repmat(n1 / norm(n1), 1, total_num);
Pgth = zeros(3 * (pairs_num + 1), total_num);
Pgth(1:3, :) = 1 ./ (n1' * q1) .* q1;
for i = 1:pairs_num
    r = r_all(i, :);
    R = rodrigues(r);
    t = t_all(i, :);
    H12 = R + t' * n1';
    q2 = H12 * q1;
    
    % double check if the depth is non-negative
    while any(q2(3, :) <= 0)
        r = randn(1, 3);
        R = rodrigues(r);
        H12 = R + t' * n1';
        q2 = H12 * q1;
    end
    r_all(i, :) = r;
    R_set{i} = R;
    q2(1, :) = q2(1, :) ./ q2(3, :);
    q2(2, :) = q2(2, :) ./ q2(3, :);
    q2(3, :) = ones(1, total_num);
    qgth{i + 1} = q2;
    I2u(i, :) = q2(1, :);
    I2v(i, :) = q2(2, :);
    
    n2 = R * n1;
    a2 = (1 + n1' * R' * t');
    H21 = R' - ((R' * t') * n2') / a2;
    Ngth((3 * i + 1):(3 * i + 3), :) = -repmat(n2 / norm(n2), 1, total_num);
    Pgth((3 * i + 1):(3 * i + 3), :) = a2 ./ (n2' * q2) .* q2;    
    [~, pu, pv, puu, puv, pvv]=homderivs(H21, q2(1:2, :));
    J21a(i, :) = pu(1, :);
    J21b(i, :) = pu(2, :);
    J21c(i, :) = pv(1, :);
    J21d(i, :) = pv(2, :);
    H21uua(i, :) = puu(1, :);
    H21uub(i, :) = puu(2, :);
    H21uva(i, :) = puv(1, :);
    H21uvb(i, :) = puv(2, :);
    H21vva(i, :) = pvv(1, :);
    H21vvb(i, :) = pvv(2, :);
    determinant = J21a(i, :) .* J21d(i, :) - J21b(i, :) .* J21c(i, :);
    J12a(i, :) = J21d(i, :) ./ determinant;
    J12b(i, :) = -J21b(i, :) ./ determinant;
    J12c(i, :) = -J21c(i, :) ./ determinant;
    J12d(i, :) = J21a(i, :) ./ determinant;
end

save('./warps_data/warps_plane_trial_critical.mat', 'H21uua', 'H21uub', 'H21uva', 'H21uvb', 'H21vva', 'H21vvb', ...
    'I1u', 'I1v', 'I2u', 'I2v', 'J21a', 'J21b', 'J21c', 'J21d', 'J12a', 'J12b', 'J12c', 'J12d', ...
    'Ngth', 'Pgth', 'qgth');
