function [na, nb, na_, nb_] = compute_normal(vec_W, vec_W_invt)
%%% n3 should be nonzero

M = cell(3);
for i = 1:3
    for j = 1:3
        M{i, j} = vec_W_invt{1, i} .* vec_W_invt{1, j} + vec_W_invt{2, i} .* vec_W_invt{2, j} + vec_W_invt{3, i} .* vec_W_invt{3, j};
    end
end

a2 = -(M{1, 1} + M{2, 2} + M{3, 3});
a1 = M{1, 1} .* M{2, 2} + M{1, 1} .* M{3, 3} + M{2, 2} .* M{3, 3} - (M{1, 2} .^ 2 + M{1, 3} .^ 2 + M{2, 3} .^ 2);
a0 = M{1, 2} .^ 2 .* M{3, 3} + M{1, 3} .^ 2 .* M{2, 2} + M{2, 3} .^ 2 .* M{1, 1} - M{1, 1} .* M{2, 2} .* M{3, 3} - 2 * M{1, 2} .* M{1, 3} .* M{2, 3};

Q = (3 * a1 - a2 .^ 2) / 9;
R = (9 * a2 .* a1 - 27 * a0 - 2 * a2 .^ 3) / 54;
S = (R + sqrt(Q .^ 3 + R .^ 2)) .^ (1 / 3);
T = (R - sqrt(Q .^ 3 + R .^ 2)) .^ (1 / 3);
lam2 = -a2 / 3 - (S + T) / 2 - sqrt(3) / 2 * (S - T) * (1i);

S = cell(3);
S{1, 1} = M{1, 1} ./ lam2 - 1; S{1, 2} = M{1, 2} ./ lam2; S{1, 3} = M{1, 3} ./ lam2;
S{2, 1} = M{2, 1} ./ lam2; S{2, 2} = M{2, 2} ./ lam2 - 1; S{2, 3} = M{2, 3} ./ lam2;
S{3, 1} = M{3, 1} ./ lam2; S{3, 2} = M{3, 2} ./ lam2; S{3, 3} = M{3, 3} ./ lam2 - 1; 

epsilon = -(S{2, 1} .* S{3, 3} - S{2, 3} .* S{3, 1});

sign_ep = sign(epsilon);
sign_ep(sign_ep == 0) = 1;

tem1 = sign_ep .* sqrt(S{1, 3} .^ 2 - S{1, 1} .* S{3, 3});
tem2 = sqrt(S{2, 3} .^ 2 - S{2, 2} .* S{3, 3});

na(:, :, 1) = S{1, 3} + tem1;
na(:, :, 2) = S{2, 3} + tem2;
na(:, :, 3) = S{3, 3};

nb(:, :, 1) = S{1, 3} - tem1;
nb(:, :, 2) = S{2, 3} - tem2;
nb(:, :, 3) = S{3, 3};

na = na ./ sqrt(na(:, :, 1) .^ 2 + na(:, :, 2) .^ 2 + na(:, :, 3) .^ 2);
nb = nb ./ sqrt(nb(:, :, 1) .^ 2 + nb(:, :, 2) .^ 2 + nb(:, :, 3) .^ 2);

na_(:, :, 1) = vec_W{1, 1} .* na(:, :, 1) + vec_W{1, 2} .* na(:, :, 2) + vec_W{1, 3} .* na(:, :, 3);
na_(:, :, 2) = vec_W{2, 1} .* na(:, :, 1) + vec_W{2, 2} .* na(:, :, 2) + vec_W{2, 3} .* na(:, :, 3);
na_(:, :, 3) = vec_W{3, 1} .* na(:, :, 1) + vec_W{3, 2} .* na(:, :, 2) + vec_W{3, 3} .* na(:, :, 3);

nb_(:, :, 1) = vec_W{1, 1} .* nb(:, :, 1) + vec_W{1, 2} .* nb(:, :, 2) + vec_W{1, 3} .* nb(:, :, 3);
nb_(:, :, 2) = vec_W{2, 1} .* nb(:, :, 1) + vec_W{2, 2} .* nb(:, :, 2) + vec_W{2, 3} .* nb(:, :, 3);
nb_(:, :, 3) = vec_W{3, 1} .* nb(:, :, 1) + vec_W{3, 2} .* nb(:, :, 2) + vec_W{3, 3} .* nb(:, :, 3);

na_ = na_ ./ sqrt(na_(:, :, 1) .^ 2 + na_(:, :, 2) .^ 2 + na_(:, :, 3) .^ 2);
nb_ = nb_ ./ sqrt(nb_(:, :, 1) .^ 2 + nb_(:, :, 2) .^ 2 + nb_(:, :, 3) .^ 2);

na = real(na);
nb = real(nb);
na_ = real(na_);
nb_ = real(nb_);

end