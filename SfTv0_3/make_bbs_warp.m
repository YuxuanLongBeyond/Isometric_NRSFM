function [ bbs, ctrlpts ] = make_bbs_warp( p, q, nC, er, KLims )

% Make a BBS warp
if nargin == 4
    umin = min(p(1,:)); umax = max(p(1,:));
    vmin = min(p(2,:)); vmax = max(p(2,:));    
elseif nargin == 5
    umin = KLims(1); umax = KLims(2);
    vmin = KLims(3); vmax = KLims(4);
end

bbs = bbs_create(umin, umax, nC, vmin, vmax, nC, size(q,1));
coloc = bbs_coloc(bbs, p(1,:), p(2,:));
lambdas = er*ones(nC-3, nC-3);
bending = bbs_bending(bbs, lambdas);

% get control points for i to j warp
cpts = (coloc'*coloc + bending) \ (coloc'*q');
ctrlpts = cpts';

end