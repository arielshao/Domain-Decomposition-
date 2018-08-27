function u=Solve1d(f,eta,a,b,gg,gd)
% SOLVE1D solves eta-Delta in 1d using finite differences
%   u=Solve1d(f,eta,a,b,gg,gd) solves the one dimensional equation
%   (eta-Delta)u=f on the domain Omega=(a,b) with Dirichlet boundary
%   conditions u=gg at x=a and u=gd at x=b using a finite
%   difference approximation with length(f) interior grid points
    
J=length(f);    
A=A1d(eta,a,b,J);         % construct 1d finite difference operator    
h=(b-a)/(J+1);
f(1)=f(1)+gg/h^2;         % add boundary conditions into rhs
f(end)=f(end)+gd/h^2;
u=A\f;
u=[gg;u;gd];              % add boundary values to solution