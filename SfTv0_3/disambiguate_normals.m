function [ N ] = disambiguate_normals( N1, N2, Nd )

N = zeros(3,size(Nd,2));
% Disambiguate the normals:
for i = 1: size(N1,2)

    n1 = N1(:,i);
    n2 = N2(:,i);
    nd = Nd(:,i);

    dot1 = n1'*nd;
    dot2 = n2'*nd;

    if abs(dot1)>abs(dot2)
        N(:,i) = -sign(n1(3))*n1;
    else
        N(:,i) = -sign(n2(3))*n2;
    end
%     N(:,i) = abs(N(:,i)).*sign(nd);
end
% N1 = Nd;