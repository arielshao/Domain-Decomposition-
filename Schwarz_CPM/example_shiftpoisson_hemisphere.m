%% Solving a shifted poisson's equation on a hemisphere

%  This example solves a shifted poisson equation on a hemisphere 
%  with dirichlet boundary condition.

%  The PDE is
%  $$-lap u+c*u = f$$ in \Omega 
%  $$       u = g$$ on \partial\Omega
%  where \Omega is a hemisphere and c is a non-negative shift.


% The method of manufactured solutions is used to construct
% an exact solution and thus determine f and g.



% This code deals with the boundary condition using
% linear extrapolation, that is,
% u(x_g)=2u(cp(x_g))-u(cpbar(x_g)) where cpbar=cp(2cp(x)-x)


% The method of manufactured solutions is used to construct
% an exact solution and thus determine f and g.

% adjust as appropriate
addpath('../cp_matrices');
addpath('../surfaces');

%% Manufactured solution
c=5;
uexact=@(th, phi)cos(3*th).*sin(phi+pi/2).^3.*(9*cos(phi+pi/2).^2-1);
f=@(th, phi, r) (30*cos(3*th).*sin(phi+pi/2).^3.*(9*cos(phi+pi/2).^2-1) )./r.^2 + c*uexact(th,phi);
g=uexact;



%% Build the mesh grid 

dx = 0.05/2;  % grid spacing
R=sqrt(2); %Radius of the circle

pad = 5;
x1d=(-R-pad*dx):dx:(R+pad*dx);
y1d=(-R-pad*dx):dx:(R+pad*dx);
z1d=(-R-pad*dx):dx:(R+pad*dx);

[xx,yy,zz]=meshgrid(x1d,y1d,z1d);
%% Banding
dim = 3;
order = 2;
p = 3;
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));

%% Closest point representations of the hemisphere
 
%closest points of the grid points
R = sqrt(2);
cpf = @(x,y,z) cpHemisphere(x,y,z,R);  
[cpx, cpy, cpz, dist,bdy] = cpf(xx, yy, zz);
  

%find the band
band = find(abs(dist) <= bw*dx);

% store closest points in the band;
cpx = cpx(band); cpy = cpy(band); cpz = cpz(band);
x = xx(band); y = yy(band); z = zz(band);
bdy=bdy(band);

[th, phi, r] = cart2sph(x,y,z);
[cpth, cpphi, cpr] = cart2sph(cpx,cpy,cpz);
  
% 3D scatter plot of the cp points
figure(1);clf;
scatter3(cpx, cpy, cpz, 40, 'filled');
hold on;
axis equal
xlabel('x'); ylabel('y'); zlabel('z');
title ('Closest point representation of the hemisphere');

%% Build the interpolating matrices and the elliptic operator
% Build the interpolating matrices
E1 = interp3_matrix(x1d, y1d, z1d, cpx, cpy, cpz, 1, band);
E = interp3_matrix(x1d, y1d, z1d, cpx, cpy, cpz, 3,band);
  
% Create Laplacian matrix
L = laplacian_3d_matrix(x1d,y1d,z1d,order,band);
 
% Penalty parameter
gamma = 2*dim/(dx^2);
I = speye(size(L));
 
% Modified elliptic operator
M=-(E1*L-gamma*(I-E))+c*I;

%% Deal with the boundary
% Compute the modified closest points which will be used to deal with boundary
[cpbarx, cpbary, cpbarz] = cpbar_3d(x, y, z, cpf);

% Build the interpolating matrix for ghost points
E_bdy=interp3_matrix(x1d,y1d,z1d, cpbarx(bdy), cpbary(bdy),cpbarz(bdy), p,band);

% Modidy the entries of the operator matrix
M(bdy,:)=(I(bdy,:)+E_bdy)/2;
%% Build RHS
rhs = f(cpth, cpphi,R);
    
% Trilinear interpolation from the grid points to the closest points of the
% ghost_points 
E_ghost = interp3_matrix(x1d, y1d, z1d, cpx(bdy), cpy(bdy), cpz(bdy), 1, band);

% Impose the boundary condition
 rhs(bdy) = E_ghost*g(th,phi); 

%  rhs(bdy) = g(th(bdy),phi(bdy));

% Solve using backslash
tic;
u = M\rhs;
toc
%% Construct an interpolation matrix for plotting on sphere
figure(2);clf;
% plotting grid on sphere, based on parameterization
[xp,yp,zp] = paramHemisphere(256,R);
xp1 = xp(:); yp1 = yp(:); zp1 = zp(:);
[th_plot, phi_plot, th_r] = cart2sph(xp1,yp1,zp1);
% Eplot is a matrix which interpolations data onto the plotting grid
Eplot = interp3_matrix(x1d, y1d, z1d, xp1, yp1, zp1, p, band);


uplot=Eplot*u;% Restrict the solution on the hemisphere
% sphplot = reshape(u_sphere, size(xp));
%surf(xp, yp, zp, sphplot);
scatter3(xp1,yp1,zp1, 40, uplot,'filled');
axis equal
xlabel('x'); ylabel('y'); zlabel('z');
title ('Solution of the poisson equation on the hemisphere');
  
%% Error analysis
 % The exact solution on the hemisphere

ue_in= uexact(th_plot,phi_plot);

% Compute the error in inf norm


err = max(abs(ue_in - uplot))
relerr=max(abs(ue_in - uplot))/max(abs(ue_in))

% c=5 the PDE is -lap u+5u=f


% dx         Error       Relative error    Elapsed time
% ------------------------------------------------------
% 0.4       2.12e-1        1.70e-1           0.019765s
% 0.2       4.25e-2        3.41e-2           0.088749s
% 0.1       1.28e-2        1.03e-2           0.476584s
% 0.05      3.76e-3        3.02e-3           3.254472s
% 0.025     8.25e-4        6.62e-4           29.063493s

% It is 2nd order accurate.

