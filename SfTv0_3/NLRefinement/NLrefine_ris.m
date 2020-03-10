function [ out ] = NLrefine_ris( bbs, ctrlpts, p, q, lbd_inext, lbd_bend)

if nargin<5
    lbd_inext = 50;
    lbd_bend = 50;
end
sz_inext = 50;
% Find initial depths
Qi = bbs_eval(bbs,ctrlpts,p(1,:)',p(2,:)',0,0);
Qz = Qi(3,:);

cam_mat = [1 0 0 0; 0 1 0 0; 0 0 1 0];
param_init = [Qz ctrlpts(1,:) ctrlpts(2,:) ctrlpts(3,:)];

% Perform non-linear refinement using [Brunet2010ACCV]

ctrlptso = ris_jacob_inext_depth_jacob(p,q,cam_mat,bbs,param_init,sz_inext,lbd_inext,lbd_bend);
% ctrlptso = ris_jacob_inext(p,q,cam_mat,bbs,ctrlpts,sz_inext,lbd_inext);
% ctrlptso = ris_jacob_inext_depth(p,q,cam_mat,bbs,param_init,sz_inext,lbd_inext);
out.phi.bbs = bbs;
out.phi.ctrlpts = ctrlptso;

end