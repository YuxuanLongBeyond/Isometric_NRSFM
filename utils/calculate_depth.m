function P_grid=calculate_depth(N_res,u,v,par, weight)
% nC=40;
nC = 20;
t = 1e-1; % 0.1
P_grid = zeros(3*size(u,1),size(u,2));
for i=1:size(u,1)
    idx = find(u(i,:)~=0) & find(v(i,:)~=0);
    umin=min(u(i,idx))-t;umax=max(u(i,idx))+t;
    vmin=min(v(i,idx))-t;vmax=max(v(i,idx))+t;
    bbsd = bbs_create(umin, umax, nC, vmin, vmax, nC, 1);
    colocd = bbs_coloc(bbsd, u(i,idx), v(i,idx));
    
    
%     u_nodes = linspace(umin, umax, nC - 3);
%     v_nodes = linspace(vmin, vmax, nC - 3);
%     [xv,yv]=meshgrid(u_nodes, v_nodes);
%     coef = par{i};
%     bbs = bbs_create(umin, umax, nC - 3, vmin, vmax, nC - 3, 1);
%     coloc = bbs_coloc(bbs, u(i, :), v(i, :));
%     bending = bbs_bending(bbs, 1e-4*ones(nC-6, nC-6));
%     cpts = (coloc'*coloc + bending) \ (coloc'*coef');
%     ctrlpts = cpts';
%     coefd =  bbs_eval(bbs,ctrlpts, xv(:),yv(:),0,0);
%     lambdas = reshape(coefd, nC - 3, nC - 3) .*ones(nC-3, nC-3);
    
    lambdas = par(i) .*ones(nC-3, nC-3);

    bendingd = bbs_bending(bbsd, lambdas);
%     [ctrlpts3Dn]=ShapeFromNormals(bbsd,colocd,bendingd,[u(i,idx);v(i,idx);ones(1,length(u(i,idx)))],N_res(3*(i-1)+1:3*(i-1)+3,idx));
    [ctrlpts3Dn]=new_SfN(bbsd,colocd,bendingd,[u(i,idx);v(i,idx);ones(1,length(u(i,idx)))],N_res(3*(i-1)+1:3*(i-1)+3,idx), weight{i});
    mu=bbs_eval(bbsd, ctrlpts3Dn,u(i,idx)',v(i,idx)',0,0);

    mu = sign(mu) .* mu;
    
    P_grid(3*(i-1)+1:3*(i-1)+3,idx) = [u(i,idx);v(i,idx);ones(1,length(u(i,idx)))].*[mu;mu;mu];
end
