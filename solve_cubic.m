function [x1, x2] = solve_cubic(eq, f1, f2, err)
n = size(eq, 2);
x2 = zeros(6, n);
x1 = x2;
for i = 1: n
    sol = roots(flipud(eq(:, i)));
    sol(abs(imag(sol)) >= err) = NaN;
    sol = [sol; NaN * ones(6 -length(sol), 1)];
    x2(:, i) = sol;

    x1(:, i) = -polyval(flipud(f2(:, i))', sol) ./ polyval(flipud(f1(:, i))', sol);
    
end


end