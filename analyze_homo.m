clear all; close all
load('data_1_100.mat');

index = find(mask2);

H_set = cell(1, length(index));
degen_metric = degen_metric(index);
for k = 1:length(index)
    i = index(k);
    H = [h11(i), h12(i), h13(i);
        h21(i), h22(i), h23(i);
        h31(i), h32(i), h33(i)];
    sigma = svd(H);
    H = (H / sigma(2));
    H_set{k} = H;
%     degen_metric(i) = norm(sigma - ones(3, 1));
end

degen_metric

error_map1(index)
error_map2(index)


H = H_set{1, 6}

na2_collect(:, index(6))
nb2_collect(:, index(6))



