function A=A1d(eta,a,b,J)
% A1D one dimensional finite difference approximation
%   A=A1d(eta,a,b,J) computes a sparse finite difference approximation
%   of the one dimensional operator eta-Delta on the domain
%   Omega=(a,b) using J interior points

h=(b-a)/(J+1); e=ones(J,1);
A=spdiags([-e/h^2 (eta+2/h^2)*e -e/h^2],[-1 0 1],J,J);