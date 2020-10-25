function P_grid = piecewise_integrate(N_res,u_all,v_all, par, cluster_num)

frame_num = size(u_all, 1);
num_p = size(u_all, 2);


for i = 1:frame_num
    coef = par(i);
    u = u_all(i, :);
    v = v_all(i, :);
    N = N_res(3 * (i - 1) + 1:(3 * i), :);
    
    X = [u; v; N]';
    
    indices = kmeans(X, cluster_num);
    for cluster = 1:cluster_num
        idx = indices == cluster;
        
        P = calculate_depth(N(:, idx), u(idx) ,v(idx), coef);
        P_grid(3 * (i - 1) + 1:(3 * i), idx) = P;
        
    end
    
end

end