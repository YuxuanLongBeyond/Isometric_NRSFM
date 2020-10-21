function nn = normal_refine(I1u, I1v, n)


nC = 20;
t = 1e-1; % 0.1
umin=min(I1u(1, :))-t;umax=max(I1u(1, :))+t;
vmin=min(I1v(1, :))-t;vmax=max(I1v(1, :))+t;
bbsd = bbs_create(umin, umax, nC, vmin, vmax, nC, 1);
colocd = bbs_coloc(bbsd, I1u(1, :), I1v(1, :));
lambdas = measure_smoothness(I1u(1, :), I1v(1, :), n) * ones(nC-3, nC-3);
bendingd = bbs_bending(bbsd, lambdas);
[ctrlpts3Dn]=ShapeFromNormals(bbsd,colocd,bendingd,[I1u(1, :);I1v(1, :);ones(1,length(I1u(1, :)))], n);
mu=bbs_eval(bbsd, ctrlpts3Dn,I1u(1, :)',I1v(1, :)',0,0);
P_grid = [I1u(1, :);I1v(1, :);ones(1,length(I1u(1, :)))].*[mu;mu;mu];


er = 1e-5;
t= 1e-3;
nC = 40;

umin=min(I1u(1, :))-t;umax=max(I1u(1, :))+t;
vmin=min(I1v(1, :))-t;vmax=max(I1v(1, :))+t;
bbs = bbs_create(umin, umax, nC, vmin, vmax, nC, 3);
coloc = bbs_coloc(bbs, I1u(1,:), I1v(1,:));
lambdas = er*ones(nC-3, nC-3);
bending = bbs_bending(bbs, lambdas);
cpts = (coloc'*coloc + bending) \ (coloc'*P_grid');
ctrlpts = cpts';

dqu = bbs_eval(bbs, ctrlpts, I1u(1, :)',I1v(1, :)',1,0);
dqv = bbs_eval(bbs, ctrlpts, I1u(1, :)',I1v(1, :)',0,1);

nu = [dqu(1,:)./sqrt(sum(dqu.^2));dqu(2,:)./sqrt(sum(dqu.^2));dqu(3,:)./sqrt(sum(dqu.^2))];
nv = [dqv(1,:)./sqrt(sum(dqv.^2));dqv(2,:)./sqrt(sum(dqv.^2));dqv(3,:)./sqrt(sum(dqv.^2))];
nn = -cross(nu,nv);

end