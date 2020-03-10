% Shape from Template example with synthetic data
%   This code is partially based on the work of  
%   [Bartoli et. al 2012] On Template-Based Reconstruction from a Single View:
%   Analytical Solutions and Proofs of Well-Posedness for  Developable,
%   Isometric and Conformal Surfaces
%   Method 'CPB's are partially based on [Chhatkuli et. al 2014]
%   Stable Template-Based Isometric 3D Reconstruction
%   in All Imaging Conditions by Linear Least-Squares.
%   Non Linear Refinement methods are based on [Brunet et. al 2014]
%   Monocular template-based 3D surface reconstruction: Convex inextensible
%   and nonconvex isometric methods
%   Method 'CPBC's are partially based on [Chhatkuli et. al 2015] submitted
%   to TPAMI
%   (c) 2013, Adrien Bartoli and Daniel Pizarro. dani.pizarro@gmail.com

% SfT is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
% 
% SfT is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
% for more details.
% 
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
clear all;
close all;
addpath('BBS');
addpath('NLRefinement');

% Script for synthetic data SfT test against number of correspondences.

load scene1s-factor1exp1N100sigma1.mat;        
p = data.p;
q = data.q;
P2 = data.P2;
indices = data.indices;

%% SfT example:
% Surface reconstruction parameters for BBS:
options.eta.er = 5e2;
options.eta.nC = 24;
options.phi.er = 8e1;
options.phi.nC = 24;
options.maxiter = 10;

options.method = 'CPBC15I';
% outd = SfTi(proi,qw,Jw,options);
outn = SfT(p,q,options);
Qn = bbs_eval(outn.phi.bbs,outn.phi.ctrlpts,p(1,:)',p(2,:)',0,0); 

% Non linear refinement with ris library
options.phi = outn.phi;
options.method = 'CPBC15I';
% outd = SfTi(proi,qw,Jw,options);
outn = SfT(p,q,options);
Qn = bbs_eval(outn.phi.bbs,outn.phi.ctrlpts,p(1,:)',p(2,:)',0,0); 

figure,
plot3(Qn(1,:),Qn(2,:),Qn(3,:),'bo');
hold on;
plot3(P2(1,:),P2(2,:),P2(3,:),'go');
axis equal;
hold off;