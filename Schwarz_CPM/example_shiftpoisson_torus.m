%% Solving a Poisson's problem on a torus 
%  This example solves the poisson's equation on a 
%  torus

%  The PDE is
%  $$-lap u+c*u = f$$ in \Omega  
%  where \Omega is a sphere and c is a non-negative shift


% The method of manufactured solutions is used to construct
% an exact solution and thus determine f and g.

% adjust as appropriate
addpath('../cp_matrices');
addpath('../surfaces');


R=1;r=0.4;
%% Manufactured solution
c=1;

% uexact=@(th, phi)sin(3*th)+cos(2*phi);
% f=@(th, phi,R,r)-(9*sin(3*th)./((R+r.*cos(phi)).^2)-2*sin(phi).*sin(2*(phi))*r./(R+r.*cos(phi)))...
%     +4*cos(2*phi)./(r.^2))+ c*uexact(th,phi);
 uexact = @(th, phi) sin(3*th) + cos(2*phi);
 f = @(th, phi)  9*sin(3*th)./(R+r*cos(phi)).^2-sin(phi)./(r*(R+r*cos(phi)))*2.*sin(2*phi)+ 4*cos(2*phi)/r^2 +c*uexact(th, phi);
%% Build the mesh grid 

dx = 0.025;  %grid size
pad = 5;
x1d=((-3-pad*dx):dx:(3+pad*dx))';
y1d=x1d;
z1d=x1d;

[xx,yy,zz]=meshgrid(x1d,y1d,z1d);


%% Banding
dim = 3;
order = 2;
p = 3;
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));

%% Closest point representations of the sphere

%closest points of the grid points
cpf = @(x,y,z) cpTorus(x,y,z,R,r,[0 0 0]);
[cpx, cpy, cpz, dist] = cpf(xx, yy, zz);
%find the band
band = find(abs(dist)<= bw*dx);
% store closest points in the band;
cpx = cpx(band); cpy = cpy(band); cpz = cpz(band);
x = xx(band); y = yy(band); z = zz(band);
[th, phi] = cart2paramTorus(x,y,z,R);
[cpth, cpphi] = cart2paramTorus(cpx,cpy,cpz,R);

% 3D scatter plot of the cp points
figure(1);clf;
scatter3(cpx, cpy, cpz, 20, 'filled');
hold on;
axis equal
xlabel('x'); ylabel('y'); zlabel('z');
title ('Closest point representation of the torus');
view(-45,45)


 
%% Build the interpolating matrices and the elliptic operator
% Build the interpolating matrices
E1 = interp3_matrix(x1d, y1d, z1d, cpx, cpy, cpz, 1, band);
E = interp3_matrix(x1d, y1d, z1d, cpx, cpy, cpz, p, band);

% Create Laplacian matrix
L = laplacian_3d_matrix(x1d,y1d,z1d, order,band);

% Penalty parameter
gamma = 2*dim/(dx^2);
I = speye(size(L));
 
% Modified elliptic operator
 M=-(E1*L-gamma*(I-E))+c*I;
 
 
%% Build RHS
rhs = f(cpth,cpphi);
% Solve using backslash
tic;
u = M\rhs;
toc
%% Construct an interpolation matrix for plotting on sphere
figure(2);clf;
% plotting grid on sphere, based on parameterization
[xp,yp,zp] = paramTorus(256,R);
xp1 = xp(:); yp1 = yp(:); zp1 = zp(:);
[th_plot, phi_plot] = cart2paramTorus(xp1,yp1,zp1,R);
% Eplot is a matrix which interpolations data onto the plotting grid
Eplot = interp3_matrix(x1d, y1d, z1d, xp1, yp1, zp1, p, band);

u_torus=Eplot*u;% Restrict the solution on the torus
 
 scatter3(xp1,yp1,zp1, 40, u_torus,'filled');
 axis equal
 xlabel('x'); ylabel('y'); zlabel('z');
 colorbar;
 view(-45,45);
 title ('Numerical solution of the shifted poisson equation on a torus');
  
%% Error analysis
% The exact solution on the sphere

ue_in= uexact(th_plot, phi_plot);

% Compute the error in inf norm
 err = max(abs(ue_in - u_torus))
 relerr=max(abs(ue_in -u_torus))/max(abs(ue_in))


%% plot of the true solution
figure(3);clf;
scatter3(xp1,yp1,zp1,40,ue_in,'filled')
axis equal
xlabel('x'); ylabel('y'); zlabel('z');
colorbar;
title ('Exact solution of the shifted poisson equation on a torus');
view(-45,45)


% The torus with outer radius R=1 and inner radius r=0.4

% dx         Error       Relative error    Elapsed time
% ------------------------------------------------------
% 0.2       1.78e-1        8.94e-2           0.083730s
% 0.1       7.99e-2        4.00e-2           0.500681s
% 0.05      1.20e-2        6.00e-3           3.574961s
% 0.025     1.70e-3        8.35e-4           61.163133s

% Order of convergence? --> 3rd order??







  



