function [ctrlpts param_opt] = ris_jacob_inext_depth_jacob(corresp_2d_tpl, corresp_2d_img, ...
                                               cam_mat, ...
                                               bbs, param_init, ...
                                               sz_inext, lbd_inext, ...
                                               lbd_bend)
%RIS_JACOB_INEXT_DEPTH_JACOB Reconstruction of 3D surface using the 
%"jacobian" inextensibility constraints
%
% SYNTAX
%  ctrlpts = ris_jacob_inext_depth_jacob(corresp_2d_tpl, corresp_2d_img, ...
%                                        cam_mat, ...
%                                        bbs, param_init, ...
%                                        sz_inext, lbd_inext, ...
%                                        lbd_bend)
%
% DESCRIPTION
%  Uses a discretization of the jacobian inextensibility constraints. The
%  optimization is simply carried with "lsqnonlin". This version is similar
%  to "ris_jacob_inext_depth" except that here the jacobian of te final
%  cost function is explicitely computed.
%  Note that in the name of the function, the first "jacob" stands for the
%  "jacobian inextensibility constraints" while the second "jacob" stands
%  for the fact that the jacobian of the cost function is computed
%  explicitely.
%  This version also includes a regularization term (based on the classical
%  bending energy).
%
% INPUT ARGUMENTS
%  - corresp_2d_tpl [2xn matrix]: matrix giving the points in the template.
%  - corresp_2d_img [2xn matrix]: matrix giving the points in the current
%     image.
%  - cam_mat [3x4 matrix]: matrix giving the projective camera matrix.
%  - bbs: b-spline structure created with BBS toolbox.
%  - param_init [n_corresp+3*n_ctrlpts vector]: initial parameters of the
%     problem. The first n_corresp values corresponds to the depth of the
%     data points. Then, there are 3*n_ctrlpts values which are the control
%     points of the B-splines (first, the n_ctrlpts values corresponding to
%     the x-coordinate of the control points, the n_ctrlpts values
%     corresponding to the y-coordinate, and so on.)
%  - sz_inext [scalar]: size of the grid used to discretize the jacobian
%     inextensibility constraints. These points are taken on a regular grid
%     of size sz_inext x sz_inext.
%  - lbd_inext [scalar>=0]: weight give to the inextensibility term in the
%     cost function.
%  - lbd_bend [scalar>=0]: weight given to the bending energy term in the
%     cost function.
%
% OUTPUT ARGUMENTS
%  - ctrlpts [3 x n_ctrlpts matrix]: control points of the refined B-spline
%     (directly usable with the BBS toolbox).
%
% (c)2010, Florent Brunet

% RIS is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
% 
% RIS is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
% for more details.
% 
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA


% Grid for discretizing the jacobian inextensibility
[u v] = meshgrid(linspace(bbs.umin,bbs.umax,sz_inext), ...
                 linspace(bbs.vmin,bbs.vmax,sz_inext));
pts_inext_tpl = [u(:)' ; v(:)'];

n_corresp = size(corresp_2d_tpl, 2);
n_ctrlpts = bbs.nptsu * bbs.nptsv;

% Square root of the "bending matrix"
lbds = lbd_bend * ones(bbs.nptsv-3, bbs.nptsu-3);
bendmat = bbs_bending(bbs, lbds);
[P D] = eig(full(bendmat));
Dpos = D;
Dpos(Dpos < 0) = 0;
sq_bendmat = P * sqrt(Dpos) / P;
sq_bendmat(abs(sq_bendmat)<=1.0e-5) = 0;
sq_bendmat = sparse(sq_bendmat);


opt = optimset('lsqnonlin');
opt = optimset(opt, 'Jacobian', 'on', 'DerivativeCheck', 'off', 'MaxIter', 40);
display(opt.MaxIter);

param_opt = lsqnonlin( ...
    @(p)cost( ...
        corresp_2d_tpl, corresp_2d_img, ...
        cam_mat, ...
        pts_inext_tpl, lbd_inext, ...
        bbs, sq_bendmat, p), ...
    param_init(:), ...
    [], [], opt);

ctrlpts_vect = param_opt(n_corresp+1:end);
ctrlpts = reshape(ctrlpts_vect, n_ctrlpts, 3)';


function [val jacob] = cost(corresp_2d_tpl, corresp_2d_img, cam_mat, ...
                            pts_inext_tpl, lbd_inext, ...
                            bbs, sq_bendmat, param)

n_corresp = size(corresp_2d_tpl, 2);
n_inext = size(pts_inext_tpl, 2);
n_ctrlpts = bbs.nptsu * bbs.nptsv;

val = zeros(3*n_corresp + 3*n_inext + 3*n_ctrlpts, 1);

depths = param(1:n_corresp);
ctrlpts_vect = param(n_corresp+1:end);
ctrlpts_mat = reshape(ctrlpts_vect, bbs.nptsu*bbs.nptsv, 3)';

%% 3D error

% 3D points with the warp
pts_3d = bbs_eval(bbs, ctrlpts_mat, corresp_2d_tpl(1,:), corresp_2d_tpl(2,:));

% Retroprojected points
sightlines = normc(cam_mat(1:3,1:3) \ [corresp_2d_img ; ones(1,n_corresp)]);
retro_pts_3d = repmat(depths',3,1) .* sightlines;

val(1:3*n_corresp) = vect(pts_3d - retro_pts_3d);


%% Jacobian inextensibility
val_bbs_du = bbs_eval(bbs, ctrlpts_mat, pts_inext_tpl(1,:), pts_inext_tpl(2,:), 1, 0);
val_bbs_dv = bbs_eval(bbs, ctrlpts_mat, pts_inext_tpl(1,:), pts_inext_tpl(2,:), 0, 1);

for i = 1:n_inext
    J = [val_bbs_du(:,i) val_bbs_dv(:,i)];
    %val(3*n_corresp+(4*(i-1)+1:4*(i-1)+4), 1) = lbd_inext * vect(J'*J - eye(2));
    tmp = lbd_inext * vect(J'*J - eye(2));
    val(3*n_corresp+(3*(i-1)+1:3*(i-1)+3), 1) = [tmp(1) ; tmp(4) ; 2*tmp(3)];
end

%% Bending energy
val(3*n_corresp+3*n_inext+(1:n_ctrlpts), 1) = ...
    sq_bendmat * param(n_corresp+(1:n_ctrlpts));
val(3*n_corresp+3*n_inext+n_ctrlpts+(1:n_ctrlpts), 1) = ...
    sq_bendmat * param(n_corresp+n_ctrlpts+(1:n_ctrlpts));
val(3*n_corresp+3*n_inext+2*n_ctrlpts+(1:n_ctrlpts), 1) = ...
    sq_bendmat * param(n_corresp+2*n_ctrlpts+(1:n_ctrlpts));

fprintf('cost=%f  (3D error=%f, inext=%f,  bend=%f)\n', sum(val.^2), ...
    sum(val(1:3*n_corresp).^2), ...
    sum(val(3*n_corresp+1:3*n_corresp+3*n_inext).^2), ...
    sum(val(3*n_corresp+3*n_inext+1:end).^2));

%% Jacobian of the cost function
if nargout >= 2
    %{
    % Trivial (but true) construction of the jacobian
    jacob = sparse(3*n_corresp+3*n_inext,n_corresp+3*n_ctrlpts);
    
    % 3D error part
    for i = 1:n_corresp
        for j = 1:3
            jacob(3*(i-1)+j,i) = - sightlines(j,i); %#ok<SPRIX>
            jacob(3*(i-1)+j,n_corresp+((j-1)*n_ctrlpts+1:j*n_ctrlpts)) ...
                = bbs_coloc(bbs, corresp_2d_tpl(1,i), corresp_2d_tpl(2,i));
        end
    end
    
    % Inextensibility part
    for i = 1:n_inext
        % Inextensibility jacobian for the current point (point #i)
        J = [val_bbs_du(:,i) val_bbs_dv(:,i)];
        w_x = bbs_coloc_deriv(bbs, pts_inext_tpl(1,i), pts_inext_tpl(2,i), 1, 0);
        w_y = bbs_coloc_deriv(bbs, pts_inext_tpl(1,i), pts_inext_tpl(2,i), 0, 1);
        
        tmp = ...
            [2*w_x*J(1,1) 2*w_x*J(2,1) 2*w_x*J(3,1) ; ...
             2*w_y*J(1,2) 2*w_y*J(2,2) 2*w_y*J(3,2) ; ...
             2*J(1,2)*w_x+2*J(1,1)*w_y ...
             2*J(2,2)*w_x+2*J(2,1)*w_y ...
             2*J(3,2)*w_x+2*J(3,1)*w_y];
         
        jacob(3*n_corresp+(3*(i-1)+1:3*(i-1)+3),n_corresp+1:end) = lbd_inext*tmp; %#ok<SPRIX>
    end
    %}
    
    nnz_bend = nnz(sq_bendmat);
    
    nz_max = 3*n_corresp*17 + 2*n_inext*3*16 + n_inext*3*32 + 3*nnz_bend;
    ii = zeros(1, nz_max);
    jj = zeros(1, nz_max);
    vv = zeros(1, nz_max);
    
    k = 1;
    
    %% 3D error
    ii(k:k+3*n_corresp-1) = 1:3*n_corresp;
    jj(k:k+3*n_corresp-1) = vect(repmat(1:n_corresp,3,1))';
    vv(k:k+3*n_corresp-1) = -sightlines(:)';
    k = k + 3*n_corresp;
    
    coloc = bbs_coloc(bbs, corresp_2d_tpl(1,:), corresp_2d_tpl(2,:));
    for i = 1:n_corresp
        tmp = coloc(i,:);
        [~, cj, cv] = find(tmp);
        
        ii(k:k+numel(cv)-1) = (3*(i-1)+1) * ones(1,numel(cv));
        jj(k:k+numel(cv)-1) = n_corresp+cj;
        vv(k:k+numel(cv)-1) = cv;
        k = k+numel(cv);
        
        ii(k:k+numel(cv)-1) = (3*(i-1)+2) * ones(1,numel(cv));
        jj(k:k+numel(cv)-1) = n_corresp+n_ctrlpts+cj;
        vv(k:k+numel(cv)-1) = cv;
        k = k+numel(cv);
        
        ii(k:k+numel(cv)-1) = (3*(i-1)+3) * ones(1,numel(cv));
        jj(k:k+numel(cv)-1) = n_corresp+2*n_ctrlpts+cj;
        vv(k:k+numel(cv)-1) = cv;
        k = k+numel(cv);
    end
    
    %% Jacobian inextensibility
    coloc_x = bbs_coloc_deriv(bbs, pts_inext_tpl(1,:), pts_inext_tpl(2,:), 1, 0);
    coloc_y = bbs_coloc_deriv(bbs, pts_inext_tpl(1,:), pts_inext_tpl(2,:), 0, 1);
    
    for i = 1:n_inext
        % Inextensibility jacobian for the current point (point #i)
        J = [val_bbs_du(:,i) val_bbs_dv(:,i)];
        w_x = 2*lbd_inext*coloc_x(i,:);
        w_y = 2*lbd_inext*coloc_y(i,:);
        
        [~, xj, xv] = find(w_x);
        [~, yj, yv] = find(w_y);
        
        ii(k:k+3*numel(xv)-1) = (3*n_corresp+3*(i-1)+1) * ones(1,3*numel(xv));
        jj(k:k+3*numel(xv)-1) = [n_corresp+xj n_corresp+n_ctrlpts+xj n_corresp+2*n_ctrlpts+xj];
        vv(k:k+3*numel(xv)-1) = [xv*J(1,1)  xv*J(2,1)  xv*J(3,1)];
        k = k+3*numel(xv);
        
        ii(k:k+3*numel(yv)-1) = (3*n_corresp+3*(i-1)+2) * ones(1,3*numel(yv));
        jj(k:k+3*numel(yv)-1) = [n_corresp+yj n_corresp+n_ctrlpts+yj n_corresp+2*n_ctrlpts+yj];
        vv(k:k+3*numel(yv)-1) = [yv*J(1,2)  yv*J(2,2)  yv*J(3,2)];
        k = k+3*numel(yv);
        
        tmp = J(1,2)*w_x+J(1,1)*w_y;
        [~, tj, tv] = find(tmp);
        ii(k:k+numel(tv)-1) = (3*n_corresp+3*(i-1)+3) * ones(1,numel(tv));
        jj(k:k+numel(tv)-1) = n_corresp+tj;
        vv(k:k+numel(tv)-1) = tv;
        k = k+numel(tv);
        
        tmp = J(2,2)*w_x+J(2,1)*w_y;
        [~, tj, tv] = find(tmp);
        ii(k:k+numel(tv)-1) = (3*n_corresp+3*(i-1)+3) * ones(1,numel(tv));
        jj(k:k+numel(tv)-1) = n_corresp+n_ctrlpts+tj;
        vv(k:k+numel(tv)-1) = tv;
        k = k+numel(tv);
        
        tmp = J(3,2)*w_x+J(3,1)*w_y;
        [~, tj, tv] = find(tmp);
        ii(k:k+numel(tv)-1) = (3*n_corresp+3*(i-1)+3) * ones(1,numel(tv));
        jj(k:k+numel(tv)-1) = n_corresp+2*n_ctrlpts+tj;
        vv(k:k+numel(tv)-1) = tv;
        k = k+numel(tv);
    end
    
    %% Bending energy
    [bi, bj, bv] = find(sq_bendmat);
    
    ii(k:k+numel(bv)-1) = 3*n_corresp+3*n_inext+bi;
    jj(k:k+numel(bv)-1) = n_corresp+bj;
    vv(k:k+numel(bv)-1) = bv;
    k = k+numel(bv);
    
    ii(k:k+numel(bv)-1) = 3*n_corresp+3*n_inext+n_ctrlpts+bi;
    jj(k:k+numel(bv)-1) = n_corresp+n_ctrlpts+bj;
    vv(k:k+numel(bv)-1) = bv;
    k = k+numel(bv);
    
    ii(k:k+numel(bv)-1) = 3*n_corresp+3*n_inext+2*n_ctrlpts+bi;
    jj(k:k+numel(bv)-1) = n_corresp+2*n_ctrlpts+bj;
    vv(k:k+numel(bv)-1) = bv;
    k = k+numel(bv);
    
    jacob = sparse(ii(1:k-1), jj(1:k-1), vv(1:k-1), ...
        3*n_corresp+3*n_inext+3*n_ctrlpts, n_corresp+3*n_ctrlpts);
end
