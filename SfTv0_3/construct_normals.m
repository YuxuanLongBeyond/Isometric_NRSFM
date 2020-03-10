function [ Nd ] = construct_normals( bbs, ctrlpts, p )
% Construct normals from the direct computation of depth
% [bbs,ctrlpts] = make_bbs_warp(proi,Q,nC,options.phi.er,options.KLims);
n1x=bbs_eval(bbs, ctrlpts, p(1,:)', p(2,:)',1,0);
n1y=bbs_eval(bbs, ctrlpts, p(1,:)', p(2,:)',0,1);
Nd=zeros(3,size(p,2));

for k=1:size(n1x,2)
    nx=n1x(:,k);
    nx=nx./norm(nx);
    ny=n1y(:,k);
    ny=ny./norm(ny);
    nn=-cross(nx,ny);
    nn=nn./norm(nn);
    Nd(:,k)=nn;
end

end