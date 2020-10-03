function [na1_collect, nb1_collect, na2_collect, nb2_collect] = normal_from_warp(I1u,I1v,I2u,I2v,J21a,J21b,J21c,J21d,J12a,J12b,J12c,J12d,H21uua,H21uub,H21uva,H21uvb,H21vva,H21vvb)


% coeff of k1       % coeff of k2        % constant term
T1_12_k1 = -J21c;   T1_12_k2 = -J21d;    T1_12_c = (J12a.*H21uva + J12c.*H21uvb);%(H21vvb./J21d)/2;
T2_12_k1 = -J21a;   T2_12_k2 = -J21b;    T2_12_c = (J12b.*H21uva + J12d.*H21uvb);%(H21uub./J21b)/2;

% k1b = -T2_12 = a*k1 + b*k2 + t1;
% k2b = -T1_12 = c*k1 + d*k2 + t2;
a = -T2_12_k1; b = -T2_12_k2; c = -T1_12_k1; d = -T1_12_k2; t1 = -T2_12_c; t2 = -T1_12_c;
pairs_num = size(a, 1);
[vec_W, vec_W_invt] = compute_W(repmat(I1u, pairs_num, 1), repmat(I1v, pairs_num, 1), I2u, I2v, a, b, c, d, t1, t2);
[na, nb, na_, nb_] = compute_normal(vec_W, vec_W_invt);

na1_collect = [na(:, :, 1); na(:, :, 2); na(:, :, 3)];
nb1_collect = [nb(:, :, 1); nb(:, :, 2); nb(:, :, 3)];

na2_collect = [na_(:, :, 1); na_(:, :, 2); na_(:, :, 3)];
nb2_collect = [nb_(:, :, 1); nb_(:, :, 2); nb_(:, :, 3)];
end