function [k3 k4 k5] = solve_k3k4k5(res,a,b,c,d,t1,t2,J12a,J12b,J12c,J12d,H21uua,H21uub,H21uva,H21uvb,H21vva,H21vvb,e,e1,u,u1,v,v1)

% find indices with no solution
idx = find(res(:,1)==0);
k1_all = [res(:,1)';a.*repmat(res(:,1)',size(e,1),1) + b.*repmat(res(:,2)',size(e,1),1) + t1];
k2_all = [res(:,2)';c.*repmat(res(:,1)',size(e,1),1) + d.*repmat(res(:,2)',size(e,1),1) + t2];


% k1 = k1_all(1,:); k1(:,idx)=[];
% k2 = k2_all(1,:); k2(:,idx)=[];
 k1 = repmat(k1_all(1,:) ,size(e,1),1);
 k2 = repmat(k2_all(1,:) ,size(e,1),1);
 k1(:,idx)=[]; k2(:,idx)=[];
k1b = k1_all(2:end,:); k1b(:,idx)=[];
k2b = k2_all(2:end,:); k2b(:,idx)=[];

% second order derivatives k3 k4 k5
Ab = (-u+ (1+u.^2).*k1b + u.*v.*k2b);    Bb = (-v + (1+v.^2).*k2b + u.*v.*k1b);
A = (-u1+ (1+u1.^2).*k1 + u1.*v1.*k2);   B = (-v1 + (1+v1.^2).*k2 + u1.*v1.*k1);


% T1b = J12a*(T1*a*a+T3*(a*b+b*a)+T5*b*b) + J12c*(T2*a*a+T4*(a*b+b*a)+T6*b*b) + J12a*Huu21a + J12c*Huu21b;
% T2b = J12b*(T1*a*a+T3*(a*b+b*a)+T5*b*b) + J12d*(T2*a*a+T4*(a*b+b*a)+T6*b*b) + J12b*Huu21a + J12d*Huu21b;
% 
% T3b = J12a*(T1*a*c+T3*(a*d+b*c)+T5*b*d) + J12c*(T2*a*c+T4*(a*d+b*c)+T6*b*d) + J12a*Huv21a + J12c*Huv21b;
% T4b = J12b*(T1*a*c+T3*(a*d+b*c)+T5*b*d) + J12d*(T2*a*c+T4*(a*d+b*c)+T6*b*d) + J12b*Huv21a + J12d*Huv21b;
% 
% T5b = J12a*(T1*c*c+T3*(c*d+c*d)+T5*d*d) + J12c*(T2*c*c+T4*(c*d+c*d)+T6*d*d) + J12a*Huu21a + J12c*Huu21b;
% T6b = J12b*(T1*c*c+T3*(c*d+c*d)+T5*d*d) + J12d*(T2*c*c+T4*(c*d+c*d)+T6*d*d) + J12b*Huu21a + J12d*Huu21b;

%second order derivatives on surface i written in terms of k3 k4 k5
% k3b = ((T1b+2*k1b)*Ab + T2b*Bb)*(1/(Ab^2+Bb^2));
% k4b = ((T3b+k2b)*Ab + (T4b+k1b)*Bb)*(1/(Ab^2+Bb^2));
% k5b = ((T6b+2*k2b)*Bb+ T5b*Ab)*(1/(Ab^2+Bb^2));

% coeff of x3
k3b100 = (Ab.*(A.*J12a.*a.^2+B.*J12c.*a.^2)+Bb.*(A.*J12b.*a.^2+B.*J12d.*a.^2))./(Ab.^2+Bb.^2);
k4b100 = (Ab.*(A.*J12a.*a.*c+B.*J12c.*a.*c)+Bb.*(A.*J12b.*a.*c+B.*J12d.*a.*c))./(Ab.^2+Bb.^2);
k5b100 = (Ab.*(A.*J12a.*c.^2+B.*J12c.*c.^2)+Bb.*(A.*J12b.*c.^2+B.*J12d.*c.^2))./(Ab.^2+Bb.^2);

% coeff of x4
k3b010 = (Ab.*(A.*J12a.*a.*b.*2.0+B.*J12c.*a.*b.*2.0)+Bb.*(A.*J12b.*a.*b.*2.0+B.*J12d.*a.*b.*2.0))./(Ab.^2+Bb.^2);
k4b010 = (Ab.*(A.*J12a.*(a.*d+b.*c)+B.*J12c.*(a.*d+b.*c))+Bb.*(A.*J12b.*(a.*d+b.*c)+B.*J12d.*(a.*d+b.*c)))./(Ab.^2+Bb.^2);
k5b010 = (Ab.*(A.*J12a.*c.*d.*2.0+B.*J12c.*c.*d.*2.0)+Bb.*(A.*J12b.*c.*d.*2.0+B.*J12d.*c.*d.*2.0))./(Ab.^2+Bb.^2);

% coeff of x5
k3b001 = (Ab.*(A.*J12a.*b.^2+B.*J12c.*b.^2)+Bb.*(A.*J12b.*b.^2+B.*J12d.*b.^2))./(Ab.^2+Bb.^2);
k4b001 = (Ab.*(A.*J12a.*b.*d+B.*J12c.*b.*d)+Bb.*(A.*J12b.*b.*d+B.*J12d.*b.*d))./(Ab.^2+Bb.^2);
k5b001 = (Ab.*(A.*J12a.*d.^2+B.*J12c.*d.^2)+Bb.*(A.*J12b.*d.^2+B.*J12d.*d.^2))./(Ab.^2+Bb.^2);

% constant term
k3b000 = -(Bb.*(J12b.*(a.^2.*k1.*2.0+a.*b.*k2.*2.0)+J12d.*(b.^2.*k2.*2.0+a.*b.*k1.*2.0)-H21uua.*J12b-H21uub.*J12d)-Ab.*(k1b.*2.0-J12a.*(a.^2.*k1.*2.0+a.*b.*k2.*2.0)-J12c.*(b.^2.*k2.*2.0+a.*b.*k1.*2.0)+H21uua.*J12a+H21uub.*J12c))./(Ab.^2+Bb.^2);
k4b000 = (Ab.*(k2b-J12a.*(k2.*(a.*d+b.*c)+a.*c.*k1.*2.0)-J12c.*(k1.*(a.*d+b.*c)+b.*d.*k2.*2.0)+H21uva.*J12a+H21uvb.*J12c)+Bb.*(k1b-J12b.*(k2.*(a.*d+b.*c)+a.*c.*k1.*2.0)-J12d.*(k1.*(a.*d+b.*c)+b.*d.*k2.*2.0)+H21uva.*J12b+H21uvb.*J12d))./(Ab.^2+Bb.^2);
k5b000 = -(Ab.*(J12a.*(c.^2.*k1.*2.0+c.*d.*k2.*2.0)+J12c.*(d.^2.*k2.*2.0+c.*d.*k1.*2.0)-H21uua.*J12a-H21uub.*J12c)-Bb.*(k2b.*2.0-J12b.*(c.^2.*k1.*2.0+c.*d.*k2.*2.0)-J12d.*(d.^2.*k2.*2.0+c.*d.*k1.*2.0)+H21uua.*J12b+H21uub.*J12d))./(Ab.^2+Bb.^2);

% metric tensor on surface 1 g111,g112,g122
g111 = e1.^2.*k1.^2 + 1 -2.*u1.*k1;
g112 = e1.^2.*k1.*k2 - u1.*k2 - v1.*k1;
g122 = e1.^2.*k2.^2 + 1 -2.*v1.*k2;

% metric tensor on surface i g11,g12,g22
g11 = e.^2.*k1b.^2 + 1 -2.*u.*k1b;
g12 = e.^2.*k1b.*k2b - u.*k2b - v.*k1b;
g22 = e.^2.*k2b.^2 + 1 -2.*v.*k2b;

%metric tensor derivatives at i
Db = (k1b.^2 + k2b.^2 + (1-u.*k1b-v.*k2b).^2);

% x11_u = k3b*D- k1b^2;
% x11_v = k4b*D - k1b*k2b;
% x21_u = k4b*D - k1b*k2b;
% x21_v = k5b*D - k2b^2;
% 
% x1_u = a*x11_u + b*x21_u;
% x1_v = c*x11_v + d*x21_v;
% x2_u = a*x11_u + b*x21_u;
% x2_v = c*x11_v + d*x21_v;
% 
% Gu = [(2*(u*a+v*b)*k1b^2  + 2*e*k1b*dx1_du           - 2*u*dx1_du        -2*a*k1b);...
%       (2*(u*a+v*b)*k1b*k2b + e*(k1b*dx2_du+dx1_du*k2b) - v*dx1_du-u*dx2_du -a*k2b-b*k1b);...
%       (2*(u*a+v*b)*k2b^2  + 2*e*k2b*dx2_du           - 2*v*dx2_du        -2*b*k2b)];
%  
% 
% Gv = [(2*(u*c+v*d)*k1b^2  + 2*e*k1b*dx1_dv           -2*u*dx1_dv        -2*c*k1b);...
%       (2*(u*c+v*d)*k1b*k2b + e*(k1b*dx2_dv+dx1_dv*k2b) -v*dx1_dv-u*dx2_dv -c*k2b-d*k1b);...
%       (2*(u*c+v*d)*k2b^2  + 2*e*k2b*dx2_dv           -2*v*dx2_dv        -2*d*k2b)];
%   

g11u100= u.*(Db.*a.*k3b100+Db.*b.*k4b100).*-2.0+e.*k1b.*(Db.*a.*k3b100+Db.*b.*k4b100).*2.0;
g11u010= u.*(Db.*a.*k3b010+Db.*b.*k4b010).*-2.0+e.*k1b.*(Db.*a.*k3b010+Db.*b.*k4b010).*2.0;
g11u001= u.*(Db.*a.*k3b001+Db.*b.*k4b001).*-2.0+e.*k1b.*(Db.*a.*k3b001+Db.*b.*k4b001).*2.0;
g11u000= a.*k1b.*(-2.0)-u.*(a.*(Db.*k3b000-k1b.^2)+b.*(Db.*k4b000-k1b.*k2b)).*2.0+k1b.^2.*(a.*u+b.*v).*2.0+e.*k1b.*(a.*(Db.*k3b000-k1b.^2)+b.*(Db.*k4b000-k1b.*k2b)).*2.0;

g12u100= -u.*(Db.*a.*k3b100+Db.*b.*k4b100)-v.*(Db.*a.*k3b100+Db.*b.*k4b100)+e.*(k1b.*(Db.*a.*k3b100+Db.*b.*k4b100)+k2b.*(Db.*a.*k3b100+Db.*b.*k4b100));
g12u010= -u.*(Db.*a.*k3b010+Db.*b.*k4b010)-v.*(Db.*a.*k3b010+Db.*b.*k4b010)+e.*(k1b.*(Db.*a.*k3b010+Db.*b.*k4b010)+k2b.*(Db.*a.*k3b010+Db.*b.*k4b010));
g12u001= -u.*(Db.*a.*k3b001+Db.*b.*k4b001)-v.*(Db.*a.*k3b001+Db.*b.*k4b001)+e.*(k1b.*(Db.*a.*k3b001+Db.*b.*k4b001)+k2b.*(Db.*a.*k3b001+Db.*b.*k4b001));
g12u000= -a.*k2b-b.*k1b-u.*(a.*(Db.*k3b000-k1b.^2)+b.*(Db.*k4b000-k1b.*k2b))-v.*(a.*(Db.*k3b000-k1b.^2)+b.*(Db.*k4b000-k1b.*k2b))+e.*(k1b.*(a.*(Db.*k3b000-k1b.^2)+b.*(Db.*k4b000-k1b.*k2b))+k2b.*(a.*(Db.*k3b000-k1b.^2)+b.*(Db.*k4b000-k1b.*k2b)))+k1b.*k2b.*(a.*u+b.*v).*2.0;

g22u100= v.*(Db.*a.*k3b100+Db.*b.*k4b100).*-2.0+e.*k2b.*(Db.*a.*k3b100+Db.*b.*k4b100).*2.0;
g22u010= v.*(Db.*a.*k3b010+Db.*b.*k4b010).*-2.0+e.*k2b.*(Db.*a.*k3b010+Db.*b.*k4b010).*2.0;
g22u001= v.*(Db.*a.*k3b001+Db.*b.*k4b001).*-2.0+e.*k2b.*(Db.*a.*k3b001+Db.*b.*k4b001).*2.0;
g22u000= b.*k2b.*(-2.0)-v.*(a.*(Db.*k3b000-k1b.^2)+b.*(Db.*k4b000-k1b.*k2b)).*2.0+k2b.^2.*(a.*u+b.*v).*2.0+e.*k2b.*(a.*(Db.*k3b000-k1b.^2)+b.*(Db.*k4b000-k1b.*k2b)).*2.0;

g11v100= u.*(Db.*c.*k4b100+Db.*d.*k5b100).*-2.0+e.*k1b.*(Db.*c.*k4b100+Db.*d.*k5b100).*2.0;
g11v010= u.*(Db.*c.*k4b010+Db.*d.*k5b010).*-2.0+e.*k1b.*(Db.*c.*k4b010+Db.*d.*k5b010).*2.0;
g11v001= u.*(Db.*c.*k4b001+Db.*d.*k5b001).*-2.0+e.*k1b.*(Db.*c.*k4b001+Db.*d.*k5b001).*2.0;
g11v000= c.*k1b.*(-2.0)-u.*(d.*(Db.*k5b000-k2b.^2)+c.*(Db.*k4b000-k1b.*k2b)).*2.0+k1b.^2.*(c.*u+d.*v).*2.0+e.*k1b.*(d.*(Db.*k5b000-k2b.^2)+c.*(Db.*k4b000-k1b.*k2b)).*2.0;

g12v100= -u.*(Db.*c.*k4b100+Db.*d.*k5b100)-v.*(Db.*c.*k4b100+Db.*d.*k5b100)+e.*(k1b.*(Db.*c.*k4b100+Db.*d.*k5b100)+k2b.*(Db.*c.*k4b100+Db.*d.*k5b100));
g12v010= -u.*(Db.*c.*k4b010+Db.*d.*k5b010)-v.*(Db.*c.*k4b010+Db.*d.*k5b010)+e.*(k1b.*(Db.*c.*k4b010+Db.*d.*k5b010)+k2b.*(Db.*c.*k4b010+Db.*d.*k5b010));
g12v001= -u.*(Db.*c.*k4b001+Db.*d.*k5b001)-v.*(Db.*c.*k4b001+Db.*d.*k5b001)+e.*(k1b.*(Db.*c.*k4b001+Db.*d.*k5b001)+k2b.*(Db.*c.*k4b001+Db.*d.*k5b001));
g12v000=-c.*k2b-d.*k1b-u.*(d.*(Db.*k5b000-k2b.^2)+c.*(Db.*k4b000-k1b.*k2b))-v.*(d.*(Db.*k5b000-k2b.^2)+c.*(Db.*k4b000-k1b.*k2b))+e.*(k1b.*(d.*(Db.*k5b000-k2b.^2)+c.*(Db.*k4b000-k1b.*k2b))+k2b.*(d.*(Db.*k5b000-k2b.^2)+c.*(Db.*k4b000-k1b.*k2b)))+k1b.*k2b.*(c.*u+d.*v).*2.0;

g22v100= v.*(Db.*c.*k4b100+Db.*d.*k5b100).*-2.0+e.*k2b.*(Db.*c.*k4b100+Db.*d.*k5b100).*2.0;
g22v010= v.*(Db.*c.*k4b010+Db.*d.*k5b010).*-2.0+e.*k2b.*(Db.*c.*k4b010+Db.*d.*k5b010).*2.0;
g22v001= v.*(Db.*c.*k4b001+Db.*d.*k5b001).*-2.0+e.*k2b.*(Db.*c.*k4b001+Db.*d.*k5b001).*2.0;
g22v000= d.*k2b.*(-2.0)-v.*(d.*(Db.*k5b000-k2b.^2)+c.*(Db.*k4b000-k1b.*k2b)).*2.0+k2b.^2.*(c.*u+d.*v).*2.0+e.*k2b.*(d.*(Db.*k5b000-k2b.^2)+c.*(Db.*k4b000-k1b.*k2b)).*2.0;


%metric tensor derivatives at i using pullback of 1
Da = (k1.^2 + k2.^2 + (1-u1.*k1-v1.*k2).^2);

% Gnew = [(g111*a*a + g112*a*b + g112*b*a + g122*b*b),...
%         (g111*a*c + g112*a*d + g112*b*c + g122*b*d);...
%         (g111*c*a + g112*c*b + g112*d*a + g122*d*b),...
%         (g111*c*c + g112*c*d + g112*d*c + g122*d*d)];

gp11 = (g111.*a.*a + g112.*a.*b + g112.*b.*a + g122.*b.*b);
gp12 = (g111.*a.*c + g112.*a.*d + g112.*b.*c + g122.*b.*d);
gp22 = (g111.*c.*c + g112.*c.*d + g112.*d.*c + g122.*d.*d);
    
% x1u = x3*D - k1^2;
% x1v = x4*D - k1*k2;
% x2u = x4*D - k1*k2;
% x2v = x5*D - k2^2;
% 
% Gu1 = (2*u1*k1^2 + 2*e1*k1*x1u -2*k1 -2*u1*x1u);
% Gu2 = (2*u1*k1*k2 + e1*(k1*x2u + x1u*x2) -v1*x1u -k2 -u1*x2u);
% Gu3 = (2*u1*k2^2 + 2*e1*k2*x2u -2*v1*x2u);
%   
% 
% Gv1 = (2*v1*k1^2 + 2*e1*k1*x1v -2*u1*x1v);
% Gv2 = (2*v1*k1*k2 + e1*(k1*x2v + x1v*k2) -v1*x1v -k1 -u1*x2v);
% Gv3 = (2*v1*k2^2 + 2*e1*k2*x2v -2*k2 -2*v1*x2v);   
% 
% gup1 = (Gu1*a*a + g111*H21uua*a + g111*a*H21uua +Gu2*a*b + g112*H21uua*b + g112*a*H21uub +Gu2*b*a + g112*H21uub*a + g112*b*H21uua +Gu3*b*b + g122*H21uub*b + g122*b*H21uub);
% gup2 = (Gu1*a*c + g111*H21uua*c + g111*a*H21uva +Gu2*a*d + g112*H21uua*d + g112*a*H21uvb +Gu2*b*c + g112*H21uub*c + g112*b*H21uva +Gu3*b*d + g122*H21uub*d + g122*b*H21uvb);
% gup3 = (Gu1*c*c + g111*H21uva*c + g111*c*H21uva +Gu2*c*d + g112*H21uva*d + g112*c*H21uvb +Gu2*d*c + g112*H21uvb*c + g112*d*H21uva +Gu3*d*d + g122*H21uvb*d + g122*d*H21uvb);
%  
% gvp1 = (Gv1*a*a + g111*H21uva*a + g111*a*H21uva +Gv2*a*b + g112*H21uva*b + g112*a*H21uvb +Gv2*b*a + g112*H21uvb*a + g112*b*H21uva +Gv3*b*b + g122*H21uvb*b + g122*b*H21uvb);
% gvp2 = (Gv1*a*c + g111*H21uva*c + g111*a*H21vva +Gv2*a*d + g112*H21uva*d + g112*a*H21vvb +Gv2*b*c + g112*H21uvb*c + g112*b*H21vva +Gv3*b*d + g122*H21uvb*d + g122*b*H21vvb);
% gvp3 = (Gv1*c*c + g111*H21vva*c + g111*c*H21vva +Gv2*c*d + g112*H21vva*d + g112*c*H21vvb +Gv2*d*c + g112*H21vvb*c + g112*d*H21vva +Gv3*d*d + g122*H21vvb*d + g122*d*H21vvb);
 
gup1100 =  -a.^2.*(Da.*u1.*2.0-Da.*e1.*k1.*2.0)-a.*b.*(Da.*v1-Da.*e1.*k2).*2.0;
gup1010 = -b.^2.*(Da.*v1.*2.0-Da.*e1.*k2.*2.0)-a.*b.*(Da.*u1-Da.*e1.*k1).*2.0;
gup1001 = zeros(size(u1));
gup1000 = -a.^2.*(k1.*2.0+e1.*k1.^3.*2.0-k1.^2.*u1.*4.0)+b.^2.*(k2.^2.*u1.*2.0+k1.*k2.*v1.*2.0-e1.*k1.*k2.^2.*2.0)+H21uua.*a.*g111.*2.0+H21uub.*a.*g112.*2.0+H21uua.*b.*g112.*2.0+H21uub.*b.*g122.*2.0-a.*b.*(k2+e1.*(k1.^2.*k2+k1.^2.*k2)-k1.^2.*v1-k1.*k2.*u1.*3.0).*2.0;


gup2100 = -a.*c.*(Da.*u1.*2.0-Da.*e1.*k1.*2.0)-a.*d.*(Da.*v1-Da.*e1.*k2)-b.*c.*(Da.*v1-Da.*e1.*k2);
gup2010 = -a.*d.*(Da.*u1-Da.*e1.*k1)-b.*c.*(Da.*u1-Da.*e1.*k1)-b.*d.*(Da.*v1.*2.0-Da.*e1.*k2.*2.0);
gup2001 = zeros(size(u1));
gup2000 = H21uva.*a.*g111+H21uvb.*a.*g112+H21uva.*b.*g112+H21uvb.*b.*g122+H21uua.*c.*g111+H21uub.*c.*g112+H21uua.*d.*g112+H21uub.*d.*g122-a.*c.*(k1.*2.0+e1.*k1.^3.*2.0-k1.^2.*u1.*4.0)-a.*d.*(k2+e1.*(k1.^2.*k2+k1.^2.*k2)-k1.^2.*v1-k1.*k2.*u1.*3.0)-b.*c.*(k2+e1.*(k1.^2.*k2+k1.^2.*k2)-k1.^2.*v1-k1.*k2.*u1.*3.0)+b.*d.*(k2.^2.*u1.*2.0+k1.*k2.*v1.*2.0-e1.*k1.*k2.^2.*2.0);


gup3100 = -c.^2.*(Da.*u1.*2.0-Da.*e1.*k1.*2.0)-c.*d.*(Da.*v1-Da.*e1.*k2).*2.0;
gup3010 =  -d.^2.*(Da.*v1.*2.0-Da.*e1.*k2.*2.0)-c.*d.*(Da.*u1-Da.*e1.*k1).*2.0;
gup3001 = zeros(size(u1));
gup3000 = -c.^2.*(k1.*2.0+e1.*k1.^3.*2.0-k1.^2.*u1.*4.0)+d.^2.*(k2.^2.*u1.*2.0+k1.*k2.*v1.*2.0-e1.*k1.*k2.^2.*2.0)+H21uva.*c.*g111.*2.0+H21uvb.*c.*g112.*2.0+H21uva.*d.*g112.*2.0+H21uvb.*d.*g122.*2.0-c.*d.*(k2+e1.*(k1.^2.*k2+k1.^2.*k2)-k1.^2.*v1-k1.*k2.*u1.*3.0).*2.0;


gvp1100 = zeros(size(u1));
gvp1010 =  -a.^2.*(Da.*u1.*2.0-Da.*e1.*k1.*2.0)-a.*b.*(Da.*v1-Da.*e1.*k2).*2.0;
gvp1001 = -b.^2.*(Da.*v1.*2.0-Da.*e1.*k2.*2.0)-a.*b.*(Da.*u1-Da.*e1.*k1).*2.0;
gvp1000 =  -b.^2.*(k2.*2.0+e1.*k2.^3.*2.0-k2.^2.*v1.*4.0)+a.^2.*(k1.^2.*v1.*2.0+k1.*k2.*u1.*2.0-e1.*k1.^2.*k2.*2.0)-a.*b.*(k1-k2.^2.*u1-k1.*k2.*v1.*3.0+e1.*k1.*k2.^2.*2.0).*2.0+H21uva.*a.*g111.*2.0+H21uvb.*a.*g112.*2.0+H21uva.*b.*g112.*2.0+H21uvb.*b.*g122.*2.0;

gvp2100 = zeros(size(u1));
gvp2010 =  -a.*c.*(Da.*u1.*2.0-Da.*e1.*k1.*2.0)-a.*d.*(Da.*v1-Da.*e1.*k2)-b.*c.*(Da.*v1-Da.*e1.*k2);
gvp2001 = -a.*d.*(Da.*u1-Da.*e1.*k1)-b.*c.*(Da.*u1-Da.*e1.*k1)-b.*d.*(Da.*v1.*2.0-Da.*e1.*k2.*2.0);
gvp2000 = -a.*d.*(k1-k2.^2.*u1-k1.*k2.*v1.*3.0+e1.*k1.*k2.^2.*2.0)-b.*c.*(k1-k2.^2.*u1-k1.*k2.*v1.*3.0+e1.*k1.*k2.^2.*2.0)+H21vva.*a.*g111+H21vvb.*a.*g112+H21vva.*b.*g112+H21vvb.*b.*g122+H21uva.*c.*g111+H21uvb.*c.*g112+H21uva.*d.*g112+H21uvb.*d.*g122-b.*d.*(k2.*2.0+e1.*k2.^3.*2.0-k2.^2.*v1.*4.0)+a.*c.*(k1.^2.*v1.*2.0+k1.*k2.*u1.*2.0-e1.*k1.^2.*k2.*2.0);


gvp3100 = zeros(size(u1));
gvp3010 =  -c.^2.*(Da.*u1.*2.0-Da.*e1.*k1.*2.0)-c.*d.*(Da.*v1-Da.*e1.*k2).*2.0;
gvp3001 =  -d.^2.*(Da.*v1.*2.0-Da.*e1.*k2.*2.0)-c.*d.*(Da.*u1-Da.*e1.*k1).*2.0;
gvp3000 = -d.^2.*(k2.*2.0+e1.*k2.^3.*2.0-k2.^2.*v1.*4.0)+c.^2.*(k1.^2.*v1.*2.0+k1.*k2.*u1.*2.0-e1.*k1.^2.*k2.*2.0)-c.*d.*(k1-k2.^2.*u1-k1.*k2.*v1.*3.0+e1.*k1.*k2.^2.*2.0).*2.0+H21vva.*c.*g111.*2.0+H21vvb.*c.*g112.*2.0+H21vva.*d.*g112.*2.0+H21vvb.*d.*g122.*2.0;


%M1 = Gnew_bar(1,1).*Gnew(2,2) - Gnew_bar(2,2)*Gnew(1,1);
L1x3 = g11u100.*gp22 + g11.*gup3100 - g22u100.*gp11 - g22.*gup1100;
L1x4 = g11u010.*gp22 + g11.*gup3010 - g22u010.*gp11 - g22.*gup1010;
L1x5 = g11u001.*gp22 + g11.*gup3001 - g22u001.*gp11 - g22.*gup1001;
L10 =  g11u000.*gp22 + g11.*gup3000 - g22u000.*gp11 - g22.*gup1000;

L2x3 = g11v100.*gp22 + g11.*gvp3100 - g22v100.*gp11 - g22.*gvp1100;
L2x4 = g11v010.*gp22 + g11.*gvp3010 - g22v010.*gp11 - g22.*gvp1010;
L2x5 = g11v001.*gp22 + g11.*gvp3001 - g22v001.*gp11 - g22.*gvp1001;
L20 =  g11v000.*gp22 + g11.*gvp3000 - g22v000.*gp11 - g22.*gvp1000;

%M2 = Gnew_bar(1,2)*Gnew(2,2) - Gnew_bar(2,2)*Gnew(1,2);
L3x3 = g12u100.*gp22 + g12.*gup3100 - g22u100.*gp12 - g22.*gup2100;
L3x4 = g12u010.*gp22 + g12.*gup3010 - g22u010.*gp12 - g22.*gup2010;
L3x5 = g12u001.*gp22 + g12.*gup3001 - g22u001.*gp12 - g22.*gup2001;
L30 =  g12u000.*gp22 + g12.*gup3000 - g22u000.*gp12 - g22.*gup2000;

L4x3 = g12v100.*gp22 + g12.*gvp3100 - g22v100.*gp12 - g22.*gvp2100;
L4x4 = g12v010.*gp22 + g12.*gvp3010 - g22v010.*gp12 - g22.*gvp2010;
L4x5 = g12v001.*gp22 + g12.*gvp3001 - g22v001.*gp12 - g22.*gvp2001;
L40 =  g12v000.*gp22 + g12.*gvp3000 - g22v000.*gp12 - g22.*gvp2000;

% remove the ill-conditioned points
idb = (Ab.^2+Bb.^2 > .01);


for i=1:size(idb,2)
    
    id = find(idb(:,i) == 0);
    L1 = [L1x3(:,i),L1x4(:,i),L1x5(:,i)]; L1c = -L10(:,i); L1(id,:)=[]; L1c(id,:) = [];
    L2 = [L2x3(:,i),L2x4(:,i),L2x5(:,i)]; L2c = -L20(:,i); L2(id,:)=[]; L2c(id,:) = [];
    L3 = [L3x3(:,i),L3x4(:,i),L3x5(:,i)]; L3c = -L30(:,i); L3(id,:)=[]; L3c(id,:) = [];
    L4 = [L4x3(:,i),L4x4(:,i),L4x5(:,i)]; L4c = -L40(:,i); L4(id,:)=[]; L4c(id,:) = [];
    
    % find second derivatives using least squares
    if ~isempty(L1)
    Y = [L1;L2;L3;L4]\[L1c;L2c;L3c;L4c];
    Y = l1decode_pd(Y, [L1;L2;L3;L4], [L1;L2;L3;L4]', [L1c;L2c;L3c;L4c]);
    k3k4k5(i,:) = [Y(1),Y(2),Y(3)];
    else
    k3k4k5(i,:) = [0,0,0];
    end
end
k3 = [k3k4k5(:,1)';k3b100.*repmat(k3k4k5(:,1)',size(idb,1),1) + ...
                       k3b010.*repmat(k3k4k5(:,2)',size(idb,1),1) + ...
                       k3b001.*repmat(k3k4k5(:,3)',size(idb,1),1) + k3b000];
k4 = [k3k4k5(:,2)';k4b100.*repmat(k3k4k5(:,1)',size(idb,1),1) + ...
                       k4b010.*repmat(k3k4k5(:,2)',size(idb,1),1) + ...
                       k4b001.*repmat(k3k4k5(:,3)',size(idb,1),1) + k4b000];
k5 = [k3k4k5(:,3)';k5b100.*repmat(k3k4k5(:,1)',size(idb,1),1) + ...
                       k5b010.*repmat(k3k4k5(:,2)',size(idb,1),1) + ...
                       k5b001.*repmat(k3k4k5(:,3)',size(idb,1),1) + k5b000];                   