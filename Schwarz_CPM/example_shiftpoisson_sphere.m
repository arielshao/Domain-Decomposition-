%% Solving a Poisson's problem on a sphere 
%  This example solves the poisson's equation on a 
%  sphere

%  The PDE is
%  $$-lap u+u = f$$ in \Omega  
%  where \Omega is a sphere


% The method of manufactured solutions is used to construct
% an exact solution and thus determine f and g.

% adjust as appropriate
addpath('../cp_matrices');
addpath('../surfaces');



%% Manufactured solution
uexact= @(th, phi) cos(3*th).*sin(phi).^3.*(9*cos(phi).^2-1);
f = @(th, phi) -29 * uexact(th, phi)
% uexact = @(x,y,z) exp(-x.*(x-1).*y.*(y-1).*z.*(z-1));
% f = @(x,y,z) - uexact(x,y,z).*( ((1-2*x).*y.*(y-1).*z.*(z-1)).^2 - 2*y.*(y-1).*z.*(z-1) + ...
%                                      ((1-2*y).*x.*(x-1).*z.*(z-1)).^2 - 2*x.*(x-1).*z.*(z-1) + ...
%                                      ((1-2*z).*x.*(x-1).*y.*(y-1)).^2 - 2*x.*(x-1).*y.*(y-1) ) ... 
%                   + uexact(x,y,z);
% uexact=@(x,y,z) 0*x +2;
% f=@(x,y,z)0*x+2;



%% Build the mesh grid 

dx = 0.1;  %grid size

% Construct the grid in the embedding space
x1d = ((-1-6*dx):dx:(1+6*dx))';
y1d = ((-1-6*dx):dx:(1+6*dx))';
z1d = ((-1-6*dx):dx:(1+6*dx))';
[xx, yy, zz] = meshgrid(x1d, y1d, z1d);


%% Banding
dim = 3;
order = 2;
p = 3;
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));

%% Closest point representations of the sphere
 
%closest points of the grid points
cpf = @(x,y,z) cpSphere(x, y, z, 1, [0 0 0]);
[cpx, cpy, cpz, dist] = cpf(xx, yy, zz);
  
%find the band
band = find(abs(dist) <= bw*dx);

% store closest points in the band;
cpx = cpx(band); cpy = cpy(band); cpz = cpz(band);
x = xx(band); y = yy(band); z = zz(band);

[th, phi, r] = cart2sph(x,y,z);
[cpth, cpphi, cpr] = cart2sph(cpx,cpy,cpz);
  
% 3D scatter plot of the cp points
figure(1);clf;
scatter3(cpx, cpy, cpz, 40, 'filled');
hold on;
axis equal
xlabel('x'); ylabel('y'); zlabel('z');
title ('Closest point representation of the sphere');


 
%% Build the interpolating matrices and the elliptic operator
% Build the interpolating matrices
E1 = interp3_matrix(x1d, y1d, z1d, cpx, cpy, cpz, 1, band);
E = interp3_matrix(x1d, y1d, z1d, cpx, cpy, cpz, p, band);
  
% Create Laplacian matrix
L = laplacian_3d_matrix(x1d,y1d,z1d,order, band);
 
% Penalty parameter
gamma = 2*dim/(dx^2) ;
I = speye(size(L)) ;
 
% Modified elliptic operator
M = -(E1*L - gamma*(I-E))+I;

 
%% Build RHS
rhs =E1* f(th, phi);
% rhs=f(cpx,cpy,cpz);

% Solve using backslash
tic;
u = M\rhs;
toc
%% Construct an interpolation matrix for plotting on sphere
figure(2);clf;
% plotting grid on sphere, based on parameterization
[xp,yp,zp] = sphere(256);
xp1 = xp(:); yp1 = yp(:); zp1 = zp(:);
[th_plot, phi_plot, th_r] = cart2sph(xp1,yp1,zp1);
% Eplot is a matrix which interpolations data onto the plotting grid
Eplot = interp3_matrix(x1d, y1d, z1d, xp1, yp1, zp1, p, band);


u_sphere=Eplot*u;% Restrict the solution on the sphere
sphplot = reshape(u_sphere, size(xp));
%surf(xp, yp, zp, sphplot);
scatter3(xp1,yp1,zp1, 40, u_sphere,'filled');
axis equal
xlabel('x'); ylabel('y'); zlabel('z');
title ('Solution of the poisson equation on the sphere');
  
%% Error analysis
 % The exact solution on the sphere
 ue_in = Eplot*uexact(th,phi); 

% Compute the error in inf norm
 err = max(abs(ue_in - u_sphere))
 relerr=max(abs(ue_in -u_sphere))/max(abs(ue_in))











  



