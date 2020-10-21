function vec_W = refine_W(vec_W_invt, N1)

vec_W = cell(3, 3);

[pairs_num, num_p] = size(vec_W_invt{1, 1});

for i = 1:pairs_num
    for j = 1:num_p
    
        
        h11 = vec_W_invt{1, 1}(i, j);
        h12 = vec_W_invt{1, 2}(i, j);
        h13 = vec_W_invt{1, 3}(i, j);

        h21 = vec_W_invt{2, 1}(i, j);
        h22 = vec_W_invt{2, 2}(i, j);
        h23 = vec_W_invt{2, 3}(i, j);

        h31 = vec_W_invt{3, 1}(i, j);
        h32 = vec_W_invt{3, 2}(i, j);
        h33 = vec_W_invt{3, 3}(i, j);
        
        % homography from first to second frame
        H_hat = [h11, h12, h13;
            h21, h22, h23;
            h31, h32, h33];
        
        [U_hat,Sigma_hat,V_hat] = svd(H_hat);
        H_hat = H_hat / Sigma_hat(2, 2);
        
        n = N1(:, j); y1 = n(1) / n(3); y2 = n(2) / n(3);
%         n_cross = [0 -n(3) n(2);
%             n(3) 0 -n(1);
%             -n(2) n(1) 0];

        S_hat = H_hat' * H_hat - eye(3);
        s_hat = [S_hat(1, 1); S_hat(1, 2); S_hat(1, 3); S_hat(2, 2); S_hat(2, 3); S_hat(3, 3)];
        
        A = [0, 0, 0, 1, -2 * y2, y2 ^ 2;
            1, 0, -2 * y1, 0, 0, y1 ^ 2;
            y2 ^ 2, -2 * y1 * y2, 0, y1 ^ 2, 0, 0];
        s = s_hat - A' * ((A * A') \ (A * s_hat));
        S = [s(1), s(2), s(3);
            s(2), s(4), s(5);
            s(3), s(5), s(6)];
        [~, Sigma_square, V] = svd(S + eye(3));
        sigma = sqrt(diag(Sigma_square));
%         H = U_hat * diag(sigma) * V';

%         W = U_hat * diag(1 ./ sigma) * V';

%         T = (H_hat * V) / sigma';
%         [U_tem, ~, V_tem] = svd(T); U = U_tem * V_tem';
%         


        M = H_hat * V * diag(sigma);
        U = M * sqrtm(M' * M) ^ (-1);
        
        
        H_ = U_hat * diag(sigma) * V';
        H = U * diag(sigma) * V';
        
        W = U * diag(1 ./ sigma) * V';
        

        vec_W{1, 1}(i, j) = W(1, 1);
        vec_W{1, 2}(i, j) = W(1, 2);
        vec_W{1, 3}(i, j) = W(1, 3);
        
        vec_W{2, 1}(i, j) = W(2, 1);
        vec_W{2, 2}(i, j) = W(2, 2);
        vec_W{2, 3}(i, j) = W(2, 3);
        
        vec_W{3, 1}(i, j) = W(3, 1);
        vec_W{3, 2}(i, j) = W(3, 2);
        vec_W{3, 3}(i, j) = W(3, 3);        
    end
end


end