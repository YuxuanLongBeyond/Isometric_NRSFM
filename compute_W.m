function [vec_W, vec_W_invt] = compute_W(I1u, I1v, I2u, I2v, a, b, c, d, t1, t2)


w11 = a + t1 .* I1u;
w12 = b + t1 .* I1v;
w21 = c + t2 .* I1u;
w22 = d + t2 .* I1v;

w13 = t1; w23 = t2;
w33 = -t1 .* I2u - t2 .* I2v + 1;
w31 = -a .* I2u - c .* I2v + I1u .* w33;
w32 = -b .* I2u - d .* I2v + I1v .* w33;

% compute cofactor of W
c11 = w22 .* w33 - w32 .* w23;
c12 = w23 .* w31 - w21 .* w33;
c13 = w21 .* w32 - w22 .* w31;
c21 = w13 .* w32 - w12 .* w33;
c22 = w11 .* w33 - w13 .* w31;
c23 = w12 .* w31 - w11 .* w32;
c31 = w12 .* w23 - w13 .* w22;
c32 = w13 .* w21 - w11 .* w23;
c33 = w11 .* w22 - w12 .* w21;

det_W = w11 .* c11 + w12 .* c12 + w13 .* c13;

vec_W = cell(3);
vec_W{1, 1} = w11; vec_W{1, 2} = w12; vec_W{1, 3} = w13; 
vec_W{2, 1} = w21; vec_W{2, 2} = w22; vec_W{2, 3} = w23; 
vec_W{3, 1} = w31; vec_W{3, 2} = w32; vec_W{3, 3} = w33; 

vec_W_invt = cell(3);
vec_W_invt{1, 1} = c11 ./ det_W; vec_W_invt{1, 2} = c12 ./ det_W; vec_W_invt{1, 3} = c13 ./ det_W; 
vec_W_invt{2, 1} = c21 ./ det_W; vec_W_invt{2, 2} = c22 ./ det_W; vec_W_invt{2, 3} = c23 ./ det_W; 
vec_W_invt{3, 1} = c31 ./ det_W; vec_W_invt{3, 2} = c32 ./ det_W; vec_W_invt{3, 3} = c33 ./ det_W; 

end