function [na, nb] = compute_normal(H)
%%% n3 should be nonzero


% sigma = svd(H);
% H = H / sigma(2);

M = H' * H;

a2 = -(M(1, 1) + M(2, 2) + M(3, 3));
a1 = M(1, 1) * M(2, 2) + M(1, 1) * M(3, 3) + M(2, 2) * M(3, 3) - (M(1, 2) ^ 2 + M(1, 3) ^ 2 + M(2, 3) ^ 2);
a0 = M(1, 2) ^ 2 * M(3, 3) + M(1, 3) ^ 2 * M(2, 2) + M(2, 3) ^ 2 * M(1, 1) - M(1, 1) * M(2, 2) * M(3, 3) - 2 * M(1, 2) * M(1, 3) * M(2, 3);

Q = (3 * a1 - a2 ^ 2) / 9;
R = (9 * a2 * a1 - 27 * a0 - 2 * a2 ^ 3) / 54;
S = (R + sqrt(Q ^ 3 + R ^ 2)) ^ (1 / 3);
T = (R - sqrt(Q ^ 3 + R ^ 2)) ^ (1 / 3);
lam2 = -a2 / 3 - (S + T) / 2 - sqrt(3) / 2 * (S - T) * (1i);

S = M / lam2 - eye(3);

% epsilon = cross(h2, h1)' * cross(h2, h3) - M(1, 3);
epsilon = -(S(2, 1) * S(3, 3) - S(2, 3) * S(3, 1));
if epsilon >= 0
    sign_ep = 1;
else
    sign_ep = -1;
end

tem1 = sign_ep * sqrt(S(1, 3) ^ 2 - S(1, 1) * S(3, 3));
tem2 = sqrt(S(2, 3) ^ 2 - S(2, 2) * S(3, 3));

na = real([S(1, 3) + tem1; S(2, 3) + tem2; S(3, 3)]);
nb = real([S(1, 3) - tem1; S(2, 3) - tem2; S(3, 3)]);

na = na / norm(na);
nb = nb / norm(nb);

end