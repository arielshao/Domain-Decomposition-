function u=Solve2d(f,eta,ai,bi,gg,gd)
% SOLVE2D solves 2d eta-Delta using a finite difference approximation
%   u=Solve2d(f,eta,a,b,gg,gd) solves the two dimensional equation
%   (eta-Delta)u=f on the domain Omega=(ai*h,bi*h)x(0,1) with
%   Dirichlet boundary conditions u=gg at x=ai*h and u=gd at x=bi*h
%   and u=0 at y=0 and y=1 using a finite difference approximation
%   with interior grid points (bi-ai) times length(gg) using the same
%   mesh size h=1/(length(gg)+1) in both x and y.

nx=bi-ai-1;    
ny=length(gg);    
h=1/(length(gg)+1);
A=A2d(eta,nx,ny);                   % Build elliptic operator 
f(1:ny,1)=f(1:ny,1)+gg/h^2;         % add boundary conditions into rhs
f(1:ny,end)=f(1:ny,end)+gd/h^2;
u=A\f(:);
u=reshape(u,ny,nx);
u=[gg u gd];                        % add boundary values to solution