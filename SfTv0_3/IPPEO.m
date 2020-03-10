function [ Q, N1, N2 ] = IPPEO( q, Jw )

% IPPE using original codes from Toby.
% Calls GFOPPE function.

N1 = zeros(3,length(q));
N2 = zeros(3,length(q));
Q = zeros(3,length(q));
for i = 1: length(q)
    v = q(1:2,i);
    
    % only for homography, change jacobian for others:
    J = [Jw(1,i) Jw(2,i); Jw(3,i) Jw(4,i)];
    P = [eye(2), -v];
    
    [R1, R2, gamma] = GFOPPE(P,J);
%     N1(:,i) = R1*[0;0;-1];
%     N2(:,i) = R2*[0;0;-1];
    N1(:,i) = R1(:,end);
    N2(:,i) = R2(:,end);
%     N1(:,i) = R2*[0;0;-1];
    Q(:,i) = [q(1,i)/gamma;q(2,i)/gamma; 1/gamma];
end
