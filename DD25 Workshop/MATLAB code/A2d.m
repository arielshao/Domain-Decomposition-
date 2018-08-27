function A=A2d(eta,Jx,Jy);
% A2D finite difference approximation of eta-Delta in 2d
%   A=A2d(eta,Jx,Jy); constructs the finite difference approximation to
%   eta-Delta on a Jx x Jy grid with spacing h=1/(Jy+1) and homogeneous
%   Dirichlet conditions all around

h=1/(Jy+1); Dxx=A1d(eta,0,h*(Jx+1),Jx); Dyy=A1d(0,0,1,Jy);
A=kron(speye(size(Dxx)),Dyy)+kron(Dxx,speye(size(Dyy)));