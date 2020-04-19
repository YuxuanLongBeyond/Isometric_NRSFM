function [I1u,I1v,I2u,I2v,J21a,J21b,J21c,J21d,J12a,J12b,J12c,J12d,H21uua,H21uub,H21uva,H21uvb,H21vva,H21vvb] = create_warps(I1u,I1v,I2u,I2v,visb,par)

% create Ground truth
er = 1e-4;
t= 1e-3;
nC = 100;
p = length(I1u(1,:));
n = length(I1u(:,1))+1;
J21a = zeros(n-1,p);
J21b = J21a; J21c = J21a; J21d = J21a; J12a = J21a; J12b = J21a; J12c = J21a; J12d = J21a;
H21uua = J21a; H21uub = J21a; H21uva = J21a; H21uvb = J21a; H21vva = J21a; H21vvb = J21a;

% find the set of possible reference images
idx = find(visb(1,:)==0);
id2 = zeros(1,length(idx));
for i = 1: length(idx)
    id = find(visb(:,idx(i))>0);
    id2(i) = id(1);
end
id = unique(id2);
id = [1,id];

Iu = [I1u(1,:);I2u]; Iv = [I1v(1,:);I2v];

% calculate Eta_21 derivatives using schwarzian warps
for j = 1:length(id)
    for i = id(j)+1:n
        idx = visb(id(j),:)==1 & visb(i,:)==1;
        q1 = [Iu(id(j),idx);Iv(id(j),idx)];  q2 = [Iu(i,idx);Iv(i,idx)];
        umin = min(q2(1,:))-t; umax = max(q2(1,:))+t; vmin = min(q2(2,:))-t; vmax = max(q2(2,:))+t;
        
        bbs = bbs_create(umin, umax, nC, vmin, vmax, nC, 2);
        coloc = bbs_coloc(bbs, q2(1,:), q2(2,:));
        lambdas = er*ones(nC-3, nC-3);
        bending = bbs_bending(bbs, lambdas);
        cpts = (coloc'*coloc + bending) \ (coloc'*q1');
        ctrlpts = cpts';
        
        [xv,yv]=meshgrid(linspace(umin,umax,100),linspace(vmin,vmax,100));
        ctrlpts = optimPanalSchwarz(bbs,ctrlpts,q2',q1',[xv(:),yv(:)],par);
        ctrlpts=ctrlpts';
        
        qw2 = bbs_eval(bbs,ctrlpts,q2(1,:)',q2(2,:)',0,0);
        error=sqrt(mean((qw2(1,:)-q1(1,:)).^2+(qw2(2,:)-q1(2,:)).^2));
        
        if j > 1
            idx = visb(1,:)==0 & visb(id(j),:)==1 & visb(i,:)==1;
            dqu = bbs_eval(bbs, ctrlpts, q2(1,idx)',q2(2,idx)',1,0); dqv = bbs_eval(bbs, ctrlpts, q2(1,idx)',q2(2,idx)',0,1);
            dquv = bbs_eval(bbs,ctrlpts,q2(1,idx)',q2(2,idx)',1,1); dquu = bbs_eval(bbs, ctrlpts, q2(1,idx)',q2(2,idx)',2,0);
            dqvv = bbs_eval(bbs, ctrlpts, q2(1,idx)',q2(2,idx)',0,2);
        else
            dqu = bbs_eval(bbs, ctrlpts, q2(1,:)',q2(2,:)',1,0); dqv = bbs_eval(bbs, ctrlpts, q2(1,:)',q2(2,:)',0,1);
            dquv = bbs_eval(bbs,ctrlpts,q2(1,:)',q2(2,:)',1,1); dquu = bbs_eval(bbs, ctrlpts, q2(1,:)',q2(2,:)',2,0);
            dqvv = bbs_eval(bbs, ctrlpts, q2(1,:)',q2(2,:)',0,2);
        end
        
        J21a(i-1,idx) = dqu(1,:); J21b(i-1,idx) = dqu(2,:); J21c(i-1,idx) = dqv(1,:); J21d(i-1,idx) = dqv(2,:);
        
        J12a(i-1,idx) = dqv(2,:)./(dqu(1,:).*dqv(2,:)-dqv(1,:).*dqu(2,:)); J12b(i-1,idx) = -dqu(2,:)./(dqu(1,:).*dqv(2,:)-dqv(1,:).*dqu(2,:));
        J12c(i-1,idx) = -dqv(1,:)./(dqu(1,:).*dqv(2,:)-dqv(1,:).*dqu(2,:)); J12d(i-1,idx) = dqu(1,:)./(dqu(1,:).*dqv(2,:)-dqv(1,:).*dqu(2,:));
        
        H21uua(i-1,idx) = dquu(1,:); H21uub(i-1,idx) = dquu(2,:); H21uva(i-1,idx) = dquv(1,:);
        H21uvb(i-1,idx) = dquv(2,:); H21vva(i-1,idx) = dqvv(1,:); H21vvb(i-1,idx) = dqvv(2,:);
        
        disp(fprintf('[ETA] Internal Rep error = %f',error));
    end
end

% arrange warps according to reference and data
idx = find(visb(1,:)==0);
for i = 1: length(idx)
    id = find(visb(1:end,idx(i))>0);
    I1u(:,idx(i)) = I2u(id(1)-1,idx(i)); I1v(:,idx(i)) = I2v(id(1)-1,idx(i)); I2u(id(1)-1,idx(i)) = 0; I2v(id(1)-1,idx(i)) = 0;
    J21a(1,idx(i)) = J21a(id(1)-1,idx(i)); J21a(id(1)-1,idx(i)) = 0;
    J21b(1,idx(i)) = J21b(id(1)-1,idx(i)); J21b(id(1)-1,idx(i)) = 0;
    J21c(1,idx(i)) = J21c(id(1)-1,idx(i)); J21c(id(1)-1,idx(i)) = 0;
    J21d(1,idx(i)) = J21d(id(1)-1,idx(i)); J21d(id(1)-1,idx(i)) = 0;
    
    J12a(1,idx(i)) = J12a(id(1)-1,idx(i)); J12a(id(1)-1,idx(i)) = 0;
    J12b(1,idx(i)) = J12b(id(1)-1,idx(i)); J12b(id(1)-1,idx(i)) = 0;
    J12c(1,idx(i)) = J12c(id(1)-1,idx(i)); J12c(id(1)-1,idx(i)) = 0;
    J12d(1,idx(i)) = J12d(id(1)-1,idx(i)); J12d(id(1)-1,idx(i)) = 0;
    
    H21uua(1,idx(i)) = H21uua(id(1)-1,idx(i)); H21uua(id(1)-1,idx(i)) = 0;
    H21uub(1,idx(i)) = H21uub(id(1)-1,idx(i)); H21uub(id(1)-1,idx(i)) = 0;
    H21uva(1,idx(i)) = H21uva(id(1)-1,idx(i)); H21uva(id(1)-1,idx(i)) = 0;
    H21uvb(1,idx(i)) = H21uvb(id(1)-1,idx(i)); H21uvb(id(1)-1,idx(i)) = 0;
    H21vva(1,idx(i)) = H21vva(id(1)-1,idx(i)); H21vva(id(1)-1,idx(i)) = 0;
    H21vvb(1,idx(i)) = H21vvb(id(1)-1,idx(i)); H21vvb(id(1)-1,idx(i)) = 0;
end