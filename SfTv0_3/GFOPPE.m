function [R1,R2,gamma] = GFOPPE(P,J)
%dd = (1/z_)*P*R1_(:,1:2) - J;
H = P*P';
G = J*J';

% H11 = H(1,1);
% H12 = H(1,2);
% H22 = H(2,2);
% 
% G11 = G(1,1);
% G12 = G(1,2);
% G22 = G(2,2);


%a_ = - H12^2 + H11*H22;
%b_ =  (2*G12*H12 - G11*H22 - G22*H11);
%c_ =  - G12^2 + G11*G22;


a = -H(1,2)^2+H(1,1)*H(2,2);
b = -H(2,2)*G(1,1) - H(1,1)*G(2,2)+2*H(1,2)*G(1,2);
c = -G(1,2)^2+G(1,1)*G(2,2);
gamma2 = ((-b+sqrt(b^2-4*a*c))/(2*a));
gamma = sqrt(gamma2);

M = gamma2*H-G;
d = [sqrt(M(1,1));sqrt(M(2,2))];
d(2) = d(2)*sign(M(1,2));

v1 = -P(1,3);
v2 = -P(2,3);

s = 1/(v1^2+v2^2+1);
Sinv = s*[v2^2+1,-v1*v2,v1;-v1*v2,v1^2+1,v2;-v1,-v2,1];

% S = P;
% S(3,:) = cross(S(1,:),S(2,:));
% 
% SInv_ = inv(S);
T1 = zeros(3,3);
T2 = zeros(3,3);

T1(1:2,:) = [J,d]/gamma;
T1(3,:) = crs(T1(1,:),T1(2,:));

T2(1:2,:) = [J,-d]/gamma;
T2(3,:) = crs(T2(1,:),T2(2,:));

R1 = Sinv*T1;
R2 = Sinv*T2;

function v3 = crs(v1,v2)
v3 = zeros(3,1);
v3(1) = v1(2)*v2(3)-v1(3)*v2(2);
v3(2) = v1(3)*v2(1)-v1(1)*v2(3);
v3(3) = v1(1)*v2(2)-v1(2)*v2(1);
