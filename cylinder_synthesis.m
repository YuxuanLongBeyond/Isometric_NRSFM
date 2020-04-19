clear all
close all

n = 10;
rmin = 2;
rmax = 11;
rot_angle = 0.5;

u = linspace(-pi/2,pi/2,20);%-1*pi/4:0.30:pi/4;
u = repmat(u,length(u),1);
aa = linspace(1,3,20);%0:2*0.095:1;
aa = repmat(aa,length(aa),1);
v = aa';
num= linspace(rmin,rmax,n);

% % Template
x = u;
y = v;

T=[x(:),y(:)];

% Make deformations
m = 0;
for i=length(num):-1:1
    r= num(i);
    x1 =  r.*sin(u/r);
    z1 =  r.*cos(u/r);
    %plot3(x1(:),y(:),z1(:),'*r');
    R{i}=rodrigues(rot_angle*[rand,rand,rand]);
    
    Tr{i} =[rand, rand, -r+3];
    P2{i}=R{i}*[x1(:), y(:),z1(:)]' + repmat(Tr{i}',1,length(x1(:)));
    dqu = R{i}*[(z1(:)./r)';(zeros(size(z1(:))))';(-x1(:)./r)'];
    dqv = R{i}*[(zeros(size(z1(:))))';(ones(size(z1(:))))';(zeros(size(z1(:))))'];
    dquu = R{i}*[(-z1(:)./(r^2))';(zeros(size(z1(:))))';(-x1(:)./(r^2))'];
    dquv = R{i}*[(zeros(size(z1(:))))';(zeros(size(z1(:))))';(zeros(size(z1(:))))'];
    dqvv = R{i}*[(zeros(size(z1(:))))';(zeros(size(z1(:))))';(zeros(size(z1(:))))'];
    jac{i} = [dqu; dqv];
    hu{i} = [dquu;dquv];
    hv{i} = [dquv;dqvv];
    for j=1:length(jac{i})
        nu = jac{i}(1:3,j);
        nu = nu./norm(nu);
        nv = jac{i}(4:6,j);
        nv = nv./norm(nv);
        normal = -cross(nu,nv);
        N_phi{i}(:,j) = normal./norm(normal);
    end
    %plot3(P2{i}(1,:),P2{i}(2,:),P2{i}(3,:),'*r');
end
 %hold off

% figure;
% hold on
% axis equal
for i=1:length(num)
    q{i}= [P2{i}(1,:)./P2{i}(3,:);P2{i}(2,:)./P2{i}(3,:)];%+0*randn(2,size(P2{i},2));
    %plot(q{i}(1,:),q{i}(2,:),'*b');
end
%hold off
er = 1e-3;
nC = 100;
t= 1e-3;
for i=1:length(num)
    umin = min(q{i}(1,:))-t; umax = max(q{i}(1,:))+t;
    vmin = min(q{i}(2,:))-t; vmax = max(q{i}(2,:))+t;
    
    bbss = bbs_create(umin, umax, nC, vmin, vmax, nC, 3);
    coloc = bbs_coloc(bbss, q{i}(1,:), q{i}(2,:));
    lambdas = er*ones(nC-3, nC-3);
    bending = bbs_bending(bbss, lambdas);
    
    % get control points for i to j warp
    cpts = (coloc'*coloc + bending) \ (coloc'*P2{i}');
    ctrlptss = cpts';
%         [xv,yv]=meshgrid(linspace(umin,umax,100),linspace(vmin,vmax,100));
%         ctrlpts = optimPanalSchwarz(bbs,ctrlpts,q{1}',T(:,1:2),[xv(:),yv(:)],1e-2);
%         ctrlpts=ctrlpts';
    q_1 = bbs_eval(bbss,ctrlptss,q{i}(1,:)',q{i}(2,:)',0,0);
    error=sqrt(mean((q_1(1,:)-P2{i}(1,:)).^2+(q_1(2,:)-P2{i}(2,:)).^2+ (q_1(3,:)-P2{i}(3,:)).^2));
    dqu = bbs_eval(bbss, ctrlptss, q{i}(1,:)',q{i}(2,:)',1,0);
    dqv = bbs_eval(bbss, ctrlptss, q{i}(1,:)',q{i}(2,:)',0,1);
    dquu = bbs_eval(bbss, ctrlptss, q{i}(1,:)',q{i}(2,:)',2,0);
    dquv = bbs_eval(bbss, ctrlptss, q{i}(1,:)',q{i}(2,:)',1,1);
    dqvv = bbs_eval(bbss, ctrlptss, q{i}(1,:)',q{i}(2,:)',0,2);
    jac1{i} = [dqu; dqv];
    hu1{i} = [dquu;dquv];
    hv1{i} = [dquv;dqvv];
    for j=1:length(jac{i})
        nu = jac1{i}(1:3,j);
        nu = nu./norm(nu);
        nv = jac1{i}(4:6,j);
        nv = nv./norm(nv);
        normal = -cross(nu,nv);
        N_phi1{i}(:,j) = normal./norm(normal);
    end
    disp([sprintf('[ETA] Internal Rep error = %f',error)]);
    
    
end


m = size(P2{1}, 2);
qgth = cell(1, n);
Pgth = zeros(3 * n, m);
Ngth = zeros(3 * n, m);
I2u = zeros(n - 1, m);
I2v = zeros(n - 1, m);
for i = 1:n
    qgth{i} = [q{i}; ones(1, m)];
    Pgth(((i - 1)*3+1):(3 * i), :) = P2{i};
    Ngth(((i - 1)*3+1):(3 * i), :) = N_phi1{i};
    if i == 1
        I1u = q{i}(1, :); I1v = q{i}(2, :);
    else
        I2u(i - 1, :) = q{i}(1, :); I2v(i - 1, :) = q{i}(2, :);
    end
end

visb = ones(n, m);
par = 2e-3;
I1u_all = repmat(I1u, n - 1, 1);
I1v_all = repmat(I1v, n - 1, 1);
[I1u,I1v,I2u,I2v,J21a,J21b,J21c,J21d,J12a,J12b,J12c,J12d,H21uua,H21uub,H21uva,H21uvb,H21vva,H21vvb] = create_warps(I1u_all,I1v_all,I2u,I2v,visb,par);

save('warps_cylinder_new2.mat', 'H21uua', 'H21uub', 'H21uva', 'H21uvb', 'H21vva', 'H21vvb', ...
    'I1u', 'I1v', 'I2u', 'I2v', 'J21a', 'J21b', 'J21c', 'J21d', 'J12a', 'J12b', 'J12c', 'J12d', ...
    'Ngth', 'Pgth', 'qgth');