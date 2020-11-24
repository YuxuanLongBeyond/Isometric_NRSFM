function [ctrpts3Dn]=new_SfN(bbs,coloc,bending,m1,n1)
coloc_du = bbs_coloc_deriv(bbs, m1(1,:)', m1(2,:)', 1, 0);
coloc_dv = bbs_coloc_deriv(bbs, m1(1,:)', m1(2,:)', 0, 1);
dot = sum(n1 .* m1);
% k1 = n1(1, :) ./ dot;
% k2 = n1(2, :) ./ dot;

A = [dot' .* coloc_du + n1(1, :)' .* coloc; dot' .* coloc_dv + n1(2, :)' .* coloc];

iter_max = 1;
w = 1;
for i = 1:iter_max
    [c, ~] = eigs((A' * (A ./ w)) + bending, 1, 'smallestabs');
    w = abs(A * c);
    w_min = median(w);
    if w_min < 1e-8
        break
    end
    w = w / w_min;
end
ctrpts3Dn =  c';
end