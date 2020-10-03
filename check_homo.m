clear all;close all;
% add libraries
addpath(genpath('BBS_NOOMP'));
addpath(genpath('tbxmanager'));
tbxmanager restorepath

addpath('SfTv0_3');
% addpath(genpath('gloptipoly3'));
% addpath(genpath('SeDuMi_1_3'));
addpath(genpath('schwarps'));
% addpath(genpath('sparseinv'));
addpath(genpath('utils'));
% addpath(genpath('l1magic'));


degen_thre = 0.15;
% degen_std = 1.5;
degen_std = 3;

error_thre = 0.15;
error_thre_after = 0.15;
frontal_thre = 0.02;

grid = 0;
grid_size = 20;
use_warp = 1;

frame1 = 1;
frame2 = 100;

im1 = imread(['kinect paper/', num2str(frame1,'%03.f'), '.bmp']);
im2 = imread(['kinect paper/', num2str(frame2,'%03.f'), '.bmp']);

dataset = 'Kinect_paper.mat';
dataset = ['./warps_data/', dataset];

load(dataset, 'qgth', 'Pgth', 'Ngth', 'K');

qgth = {qgth{frame1}, qgth{frame2}};
Ngth = [Ngth((3 *(frame1 - 1) + 1): (3 * frame1), :); Ngth((3 *(frame2 - 1) + 1): (3 * frame2), :)];
Pgth = [Pgth((3 *(frame1 - 1) + 1): (3 * frame1), :); Pgth((3 *(frame2 - 1) + 1): (3 * frame2), :)];


%3D Ground truth points and normalized image points
n = length(qgth);
pairs_num = n - 1;
for i=1:n
    q_n(2*(i-1)+1:2*(i-1)+2,:) = qgth{i}(1:2, :);
end


visb = ones(2, size(q_n, 2));

if grid %max(sum(visb)) == num
    % make a grid
    [I1u,I1v,I2u,I2v,visb] = create_grid(q_n,visb,grid_size);
else
    % point-wise
    I1u = repmat(q_n(1,:),n-1,1);
    I1v = repmat(q_n(2,:),n-1,1);
    I2u = q_n(3:2:2*n,:);
    I2v = q_n(4:2:2*n,:);
end


if grid
    er = 1e-5;
    t= 1e-3;
    nC = 40;
    
    umin = min([q_n(1,:), I1u])-t; umax = max([q_n(1,:), I1u])+t;
    vmin = min([q_n(2,:), I1v])-t; vmax = max([q_n(2,:), I1v])+t;
    
    bbs = bbs_create(umin, umax, nC, vmin, vmax, nC, 3);
    coloc = bbs_coloc(bbs, q_n(1,:), q_n(2,:));
    lambdas = er*ones(nC-3, nC-3);
    bending = bbs_bending(bbs, lambdas);
    cpts = (coloc'*coloc + bending) \ (coloc'*Pgth(1:3, :)');
    ctrlpts = cpts';
    
    dqu = bbs_eval(bbs, ctrlpts, I1u',I1v',1,0);
    dqv = bbs_eval(bbs, ctrlpts, I1u',I1v',0,1);
    
    nu = [dqu(1,:)./sqrt(sum(dqu.^2));dqu(2,:)./sqrt(sum(dqu.^2));dqu(3,:)./sqrt(sum(dqu.^2))];
    nv = [dqv(1,:)./sqrt(sum(dqv.^2));dqv(2,:)./sqrt(sum(dqv.^2));dqv(3,:)./sqrt(sum(dqv.^2))];
    nn = -cross(nu,nv);
    n1gth = [nn(1,:)./sqrt(sum(nn.^2));nn(2,:)./sqrt(sum(nn.^2));nn(3,:)./sqrt(sum(nn.^2))];
    p1gth = bbs_eval(bbs, ctrlpts, I1u',I1v',0,0);
    
    umin = min([q_n(3,:), I2u])-t; umax = max([q_n(3,:), I2u])+t;
    vmin = min([q_n(4,:), I2v])-t; vmax = max([q_n(4,:), I2v])+t;
    
    bbs = bbs_create(umin, umax, nC, vmin, vmax, nC, 3);
    coloc = bbs_coloc(bbs, q_n(3,:), q_n(4,:));
    lambdas = er*ones(nC-3, nC-3);
    bending = bbs_bending(bbs, lambdas);
    cpts = (coloc'*coloc + bending) \ (coloc'*Pgth(4:6, :)');
    ctrlpts = cpts';
    
    dqu = bbs_eval(bbs, ctrlpts, I2u',I2v',1,0);
    dqv = bbs_eval(bbs, ctrlpts, I2u',I2v',0,1);
    
    nu = [dqu(1,:)./sqrt(sum(dqu.^2));dqu(2,:)./sqrt(sum(dqu.^2));dqu(3,:)./sqrt(sum(dqu.^2))];
    nv = [dqv(1,:)./sqrt(sum(dqv.^2));dqv(2,:)./sqrt(sum(dqv.^2));dqv(3,:)./sqrt(sum(dqv.^2))];
    nn = -cross(nu,nv);
    n2gth = [nn(1,:)./sqrt(sum(nn.^2));nn(2,:)./sqrt(sum(nn.^2));nn(3,:)./sqrt(sum(nn.^2))];
    p2gth = bbs_eval(bbs, ctrlpts, I2u',I2v',0,0);
else
    n1gth = Ngth(1:3, :);
    n2gth = Ngth(4:6, :);
    p1gth = Pgth(1:3, :);
    p2gth = Pgth(4:6, :);
end

% 
if use_warp
    load(['kinect_schwarps_', int2str(frame1), '_', int2str(frame2), '.mat'], 'I1u','I1v','I2u','I2v','J21a','J21b','J21c','J21d','J12a','J12b','J12c','J12d','H21uua','H21uub','H21uva','H21uvb','H21vva','H21vvb')
%     load(dataset, 'I1u','I1v','I2u','I2v','J21a','J21b','J21c','J21d','J12a','J12b','J12c','J12d','H21uua','H21uub','H21uva','H21uvb','H21vva','H21vvb');
%     I1u = I1u(1, :); I1v = I1v(1, :);
%     I2u = I2u(1, :); I2v = I2v(1, :);
%     J21a = J21a(frame2 - 1, :); J21b = J21b(frame2 - 1, :);
%     J21c = J21c(frame2 - 1, :); J21d = J21d(frame2 - 1, :);
%     J12a = J12a(frame2 - 1, :); J12b = J12b(frame2 - 1, :);
%     J12c = J12c(frame2 - 1, :); J12d = J12d(frame2 - 1, :);
%     H21uua = H21uua(frame2 - 1, :); H21uub = H21uub(frame2 - 1, :);
%     H21uva = H21uva(frame2 - 1, :); H21uvb = H21uvb(frame2 - 1, :);
%     H21vva = H21vva(frame2 - 1, :); H21vvb = H21vvb(frame2 - 1, :);
else
    par = 2e-3;
    visb = ones(2, size(I1u, 2));
    [I1u,I1v,I2u,I2v,J21a,J21b,J21c,J21d,J12a,J12b,J12c,J12d,H21uua,H21uub,H21uva,H21uvb,H21vva,H21vvb] = create_warps(I1u,I1v,I2u,I2v,visb,par,1);
end
image_coord1 = K * [I1u; I1v; ones(1, length(I1u))];
image_coord1 = image_coord1 ./ image_coord1(3, :);
image_coord2 = K * [I2u; I2v; ones(1, length(I2u))];
image_coord2 = image_coord2 ./ image_coord2(3, :);

% coeff of k1       % coeff of k2        % constant term
T1_12_k1 = -J21c;   T1_12_k2 = -J21d;    T1_12_c = (J12a.*H21uva + J12c.*H21uvb);%(H21vvb./J21d)/2;
T2_12_k1 = -J21a;   T2_12_k2 = -J21b;    T2_12_c = (J12b.*H21uva + J12d.*H21uvb);%(H21uub./J21b)/2;

% k1b = -T2_12 = a*k1 + b*k2 + t1;
% k2b = -T1_12 = c*k1 + d*k2 + t2;
a = -T2_12_k1; b = -T2_12_k2; c = -T1_12_k1; d = -T1_12_k2; t1 = -T2_12_c; t2 = -T1_12_c;

[vec_W, vec_W_invt] = compute_W(repmat(I1u, pairs_num, 1), repmat(I1v, pairs_num, 1), I2u, I2v, a, b, c, d, t1, t2);
[na, nb, na_, nb_, sigma] = compute_normal(vec_W, vec_W_invt);

sigma = [sigma(:, :, 1); sigma(:, :, 2); sigma(:, :, 3)];
degen_metric = (sqrt(sum((sigma - 1) .^ 2)) + sqrt(sum((1 ./ sigma - 1) .^ 2))) / 2;
% degen_tem = degen_metric;
% h11 = vec_W{1, 1}; h12 = vec_W{2, 1}; h13 = vec_W{3, 1};
% h21 = vec_W{1, 2}; h22 = vec_W{2, 2}; h23 = vec_W{3, 2};
% h31 = vec_W{1, 3}; h32 = vec_W{2, 3}; h33 = vec_W{3, 3};
% degen_metric = zeros(1, size(I1u, 2));
% 
% for i = 1:size(I1u, 2)
%     H = [h11(i), h12(i), h13(i);
%         h21(i), h22(i), h23(i);
%         h31(i), h32(i), h33(i)];
%     s = svd(H);
%     s = (s / s(2));
%     degen_metric(i) = norm(s - ones(3, 1));
% end
mask_degen = (degen_metric < degen_thre) & (degen_metric < (mean(degen_metric)- degen_std * std(degen_metric)));

degen_idx = find(mask_degen);

na1_collect = [na(:, :, 1); na(:, :, 2); na(:, :, 3)];
nb1_collect = [nb(:, :, 1); nb(:, :, 2); nb(:, :, 3)];

na2_collect = [na_(:, :, 1); na_(:, :, 2); na_(:, :, 3)];
nb2_collect = [nb_(:, :, 1); nb_(:, :, 2); nb_(:, :, 3)];

[error_map1, index1] = min([sqrt(sum(cross(na1_collect, n1gth) .^ 2)); sqrt(sum(cross(nb1_collect, n1gth) .^ 2))]);
[error_map2, index2] = min([sqrt(sum(cross(na2_collect, n2gth) .^ 2)); sqrt(sum(cross(nb2_collect, n2gth) .^ 2))]);

disp('average minimum shape error of first frame, before integration')
mean(error_map1)
disp('average minimum shape error of second frame, before integration')
mean(error_map2)

for i = 1:length(index1)
    if index1(i) == 1
        N_res(1:3, i) = na1_collect(:, i);
    else
        N_res(1:3, i) = nb1_collect(:, i);
    end
    
%     % test normal integration
%     N_res(1:3, i) = n1gth(:, i);
    
    if index2(i) == 1
        N_res(4:6, i) = na2_collect(:, i);
    else
        N_res(4:6, i) = nb2_collect(:, i);
    end    
end

N_res(1:3, :) = sign(N_res(3, :)) .* N_res(1:3, :);
N_res(4:6, :) = sign(N_res(6, :)) .* N_res(4:6, :);

figure
draw_surface(p1gth, 'g')
hold on
quiver3(p1gth(1, :), p1gth(2, :), p1gth(3, :), N_res(1, :), N_res(2, :), N_res(3, :), 'r')
% quiver3(Pgth(1, :), Pgth(2, :), Pgth(3, :), Ngth(1, :) .* sign(Ngth(3, :)), Ngth(2, :) .* sign(Ngth(3, :)), Ngth(3, :) .* sign(Ngth(3, :)), 'r')
axis equal
hold off
figure
draw_surface(p2gth, 'g')
hold on
quiver3(p2gth(1, :), p2gth(2, :), p2gth(3, :), N_res(4, :), N_res(5, :), N_res(6, :), 'r')
% quiver3(Pgth(4, :), Pgth(5, :), Pgth(6, :), Ngth(4, :) .* sign(Ngth(6, :)), Ngth(5, :) .* sign(Ngth(6, :)), Ngth(6, :) .* sign(Ngth(6, :)), 'r')
axis equal

u_all = [I1u; I2u];
v_all = [I1v; I2v];

u_all(:, degen_idx) = 0;
v_all(:, degen_idx) = 0;

P_grid = calculate_depth(N_res,u_all,v_all,1e0);

% compare with ground truth
[P2,err_p] = compare_with_Pgth(P_grid,u_all,v_all,q_n,Pgth);
[N,err_n] = compare_with_Ngth(P2,q_n,Ngth);

figure();     

draw_surface(Pgth(1:3, :), 'g')
hold on
draw_surface(P2(1:3, :), 'r')
hold off
figure();
draw_surface(Pgth(4:6, :), 'g')
hold on
draw_surface(P2(4:6, :), 'r')
hold off

error_map1_after = sqrt(sum(cross(Ngth(1:3, :), N(1:3, :)) .^ 2));
error_map2_after = sqrt(sum(cross(Ngth(4:6, :), N(4:6, :)) .^ 2));

disp('average minimum shape error of first frame, after integration')
mean(error_map1_after)
disp('average minimum shape error of second frame, after integration')
mean(error_map2_after)

mask1 = error_map1 > error_thre;
mask2 = error_map2 > error_thre;

figure
hold off
subplot(2,2, 1)
imshow(im1)
hold on
plot(image_coord1(1, mask1), image_coord1(2, mask1), '+r')
% plot(image_coord1(1, mask1_after), image_coord1(2, mask1_after), '+g')
title(['normal error (reference) greater than ', num2str(error_thre)])
subplot(2,2, 2)
imshow(im2)
hold on
plot(image_coord2(1, mask2), image_coord2(2, mask2), '+r')
% plot(image_coord2(1, mask2_after), image_coord2(2, mask2_after), '+g')
title(['normal error greater than ', num2str(error_thre)])



hold off
subplot(2,2, 3)
imshow(im1)
hold on
plot(image_coord1(1, mask_degen), image_coord1(2, mask_degen), '*r')
title(['non-degeneracy (reference) less than ', num2str(degen_thre)])
subplot(2,2, 4)
imshow(im2)
hold on
plot(image_coord2(1, mask_degen), image_coord2(2, mask_degen), '*r')
title(['non-degeneracy less than ', num2str(degen_thre)])


measure_frontal1 = (N_res(1, :) .^ 2 + N_res(2, :) .^ 2) ./ N_res(3, :) .^ 2;
measure_frontal2 = (N_res(4, :) .^ 2 + N_res(5, :) .^ 2) ./ N_res(6, :) .^ 2;

% measure_frontal1 = (Ngth(1, :) .^ 2 + Ngth(2, :) .^ 2) ./ Ngth(3, :) .^ 2;
% measure_frontal2 = (Ngth(4, :) .^ 2 + Ngth(5, :) .^ 2) ./ Ngth(6, :) .^ 2;
% 

mask1_after = error_map1_after > error_thre_after;
mask2_after = error_map2_after > error_thre_after;

image_raw_coord1 = K * qgth{1};
image_raw_coord2 = K * qgth{2};
image_raw_coord1 = image_raw_coord1 ./ image_raw_coord1(3, :);
image_raw_coord2 = image_raw_coord2 ./ image_raw_coord2(3, :);

mask1_frontal = measure_frontal1 < frontal_thre;
mask2_frontal = measure_frontal2 < frontal_thre;
figure
hold off
subplot(2,2, 1)
imshow(im1)
hold on
plot(image_raw_coord1(1, mask1_after), image_raw_coord1(2, mask1_after), '*r')
title(['normal error (after integration) greater than', num2str(error_thre_after)])
subplot(2,2, 2)
imshow(im2)
hold on
plot(image_raw_coord2(1, mask2_after), image_raw_coord2(2, mask2_after), '*r')
title(['normal error (after integration) greater than', num2str(error_thre_after)])

subplot(2,2, 3)
imshow(im1)
hold on
plot(image_coord1(1, mask1_frontal), image_coord1(2, mask1_frontal), '*r')
title(['frontal parallelism (reference) less than ', num2str(frontal_thre)])
subplot(2,2, 4)
imshow(im2)
hold on
plot(image_coord2(1, mask2_frontal), image_coord2(2, mask2_frontal), '*r')
title(['frontal parallelism less than ', num2str(frontal_thre)])


% load('schwarp_1_40_error.mat');
% yd = reshape(y(1:(2 * size(q_n, 2))), 2, size(q_n, 2));
% interp_err = sum(yd .^ 2);
% 
% mask_interp = interp_err > 0.00001;
% 
% figure
% hold off
% subplot(1,2, 1)
% imshow(im1)
% hold on
% plot(image_raw_coord1(1, mask_interp), image_raw_coord1(2, mask_interp), '*r')
% title(['normal error (after integration) greater than', num2str(0.0001)])
% subplot(1,2, 2)
% imshow(im2)
% hold on
% plot(image_raw_coord2(1, mask_interp), image_raw_coord2(2, mask_interp), '*r')
% title(['normal error (after integration) greater than', num2str(0.0001)])



% num = length(index);
% dub_H = cell(1, num);
% for i = 1:num
%     id = index(i);
%     dub_H{i} = [h11(frame - 1, id), h12(frame - 1, id), h13(frame - 1, id);
%         h21(frame - 1, id), h22(frame - 1, id), h23(frame - 1, id);
%         h31(frame - 1, id), h32(frame - 1, id), h33(frame - 1, id)];
% end
% 
% for i = 1:num
%     H = dub_H{i};
%     id = index(i);
% %     x = [I1u(frame - 1, id); ]
%     x = 0; y = 0;
%     [R, t, n, d] = decompose(H,x,y);
%     
%     n = [n{1}, n{2}, n{3}, n{4}];
%     msk = find(n(3, :) > 0);
%     
%     n = n(:, msk);
%     na2(:, i) = n(:, 1); nb2(:, i) = n(:, 2);
%     na1(:, i) = R{msk(1)} * n(:, 1); nb1(:, i) = R{msk(2)} * n(:, 2);
% end
% 
% load('normals_1_40.mat');
% N_selected1 = N_res(1:3, index);
% N_not_selected1 = N_res_another(1:3, index);
% 
% N_selected2 = N_res(4:6, index);
% N_not_selected2 = N_res_another(4:6, index);
% 
% Ng1_picked = Ng1(:, index);
% Ng2_picked = Ng2(:, index);
% 
% error1 = sqrt(sum(cross(Ng1_picked, N_selected1) .^ 2));
% error2 = sqrt(sum(cross(Ng2_picked, N_selected2) .^ 2));
% 
% load('error_maps_1_40.mat')
% 
% 
% norm(N_selected1 - na1)
% norm(N_selected2 - na2)
% 
% load('warps_data/Kinect_paper.mat')
% % coeff of k1       % coeff of k2        % constant term
% T1_12_k1 = -J21c;   T1_12_k2 = -J21d;    T1_12_c = (J12a.*H21uva + J12c.*H21uvb);%(H21vvb./J21d)/2;
% T2_12_k1 = -J21a;   T2_12_k2 = -J21b;    T2_12_c = (J12b.*H21uva + J12d.*H21uvb);%(H21uub./J21b)/2;
% 
% % k1b = -T2_12 = a*k1 + b*k2 + t1;
% % k2b = -T1_12 = c*k1 + d*k2 + t2;
% a = -T2_12_k1; b = -T2_12_k2; c = -T1_12_k1; d = -T1_12_k2; t1 = -T2_12_c; t2 = -T1_12_c;
% 
% a_picked = a(frame - 1, index);
% b_picked = b(frame - 1, index);
% c_picked = c(frame - 1, index);
% d_picked = d(frame - 1, index);
% t1_picked = t1(frame - 1, index);
% t2_picked = t2(frame - 1, index);
% 
% J12a_picked = J12a(frame - 1, index);
% J12b_picked = J12b(frame - 1, index);
% J12c_picked = J12c(frame - 1, index);
% J12d_picked = J12d(frame - 1, index);
% 
% H21uva_picked = H21uva(frame - 1, index);
% H21uvb_picked = H21uvb(frame - 1, index);
% H21uva_picked = reshape(H21uva_picked, 3, 3);
% H21uvb_picked = reshape(H21uvb_picked, 3, 3);
% 
% 
% a_picked = reshape(a_picked, 3, 3);
% b_picked = reshape(b_picked, 3, 3);
% c_picked = reshape(c_picked, 3, 3);
% d_picked = reshape(d_picked, 3, 3);
% t1_picked = reshape(t1_picked, 3, 3);
% t2_picked = reshape(t2_picked, 3, 3);
% 
% J12a_picked = reshape(J12a_picked, 3, 3);
% J12b_picked = reshape(J12b_picked, 3, 3);
% J12c_picked = reshape(J12c_picked, 3, 3);
% J12d_picked = reshape(J12d_picked, 3, 3);
% 
% error1 = reshape(error1, 3,3);
% error2 = reshape(error2, 3,3);