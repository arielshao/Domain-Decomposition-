function u=Solve2dR(f,eta,ai,bi,gg,gd,p1,p2)
% SOLVE2DR solves 2d eta-Delta using Robin conditions
%   u=Solve2dR(f,eta,a,b,gg,gd,n,p1,p2) solves the two dimensional
%   equation (eta-Delta)u=f on the domain Omega=(ai*h,bi*h)x(0,1) with
%   Robin boundary conditions (dn+p1)u=gg at x=ai*h and (dn+p2)u=gd at
%   x=bi*h and u=0 at y=0 and y=1 using a finite difference
%   approximation with interior grid points (bi-ai) times length(gg)
%   using the same mesh size h=1/(length(gg)+1) in both x and y.

nx=bi-ai+1;    
ny=length(gg);    
h=1/(length(gg)+1);
A=A2d(eta,nx,ny);
A(1:ny,1:ny)=A(1:ny,1:ny)/2+p1/h*speye(ny);
A(end-ny+1:end,end-ny+1:end)=A(end-ny+1:end,end-ny+1:end)/2+p2/h*speye(ny);
f(1:ny,1)=f(1:ny,1)/2+gg/h;         % add boundary conditions into rhs
f(1:ny,end)=f(1:ny,end)/2+gd/h;
u=A\f(:);
u=reshape(u,ny,nx);