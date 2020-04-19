function Ngth = create_gth_normals(P2_n,q_n,n)

% create Ground truth
er = 1e-5;
t= 1e-3;
nC = 40;
Ngth = zeros(size(P2_n));
for i=1:n
    idx = find(q_n(2*(i-1)+1,:)~=0);
    umin = min(q_n(2*(i-1)+1,idx))-t; umax = max(q_n(2*(i-1)+1,idx))+t;
    vmin = min(q_n(2*(i-1)+2,idx))-t; vmax = max(q_n(2*(i-1)+2,idx))+t;
    
    bbs = bbs_create(umin, umax, nC, vmin, vmax, nC, 3);
    coloc = bbs_coloc(bbs, q_n(2*(i-1)+1,idx), q_n(2*(i-1)+2,idx));
    lambdas = er*ones(nC-3, nC-3);
    bending = bbs_bending(bbs, lambdas);
    cpts = (coloc'*coloc + bending) \ (coloc'*P2_n(3*(i-1)+1:3*(i-1)+3,idx)');
    ctrlpts = cpts';
    
    dqu = bbs_eval(bbs, ctrlpts, q_n(2*(i-1)+1,idx)',q_n(2*(i-1)+2,idx)',1,0);
    dqv = bbs_eval(bbs, ctrlpts, q_n(2*(i-1)+1,idx)',q_n(2*(i-1)+2,idx)',0,1);
    
    nu = [dqu(1,:)./sqrt(sum(dqu.^2));dqu(2,:)./sqrt(sum(dqu.^2));dqu(3,:)./sqrt(sum(dqu.^2))];
    nv = [dqv(1,:)./sqrt(sum(dqv.^2));dqv(2,:)./sqrt(sum(dqv.^2));dqv(3,:)./sqrt(sum(dqv.^2))];
    nn = -cross(nu,nv);
    Ngth(3*(i-1)+1:3*(i-1)+3,idx) = [nn(1,:)./sqrt(sum(nn.^2));nn(2,:)./sqrt(sum(nn.^2));nn(3,:)./sqrt(sum(nn.^2))];
end

