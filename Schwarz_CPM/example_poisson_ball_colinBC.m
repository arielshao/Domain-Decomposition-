%% Solving a Poisson's problem on a solid ball with Dirichlet B.C
%  This example solves the poisson's equation on a solid ball 
%  with dirichlet boundary condition.

%  The PDE is
%  $$lap u = f$$ in \Omega  
%  $$u=g$$       on\partial\Omega 
%  where \Omega is the solid ball


% The method of manufactured solutions is used to construct
% an exact solution and thus determine f and g.


% This code works well for homogenous b.c.

% adjust as appropriate
addpath('../cp_matrices');
addpath('../surfaces');



 %% Manufactured solution
%   uexact = @(x,y,z) (exp(sin(4*x)) + x.*sin(5*y))/2 + z;
%   uxexact = @(x,y,z) 2 * exp(sin (4 * x)) .* cos(4 * x) + sin(5 * y) / 2;
%   uyexact = @(x,y,z) 5 * x .* cos(5 * y) / 2;
%   uzexact = @(x,y,z) ones(size(x));
%   f = @(x,y,z) -25 * x .* sin (5 * y) / 2 + 8 * (-sin (4 * x) + cos (4 * x) .^ 2) .* exp (sin (4 * x));
% %   g=uexact;

uexact=@(x,y,z) (1/6)* (x.^2+y.^2+z.^2-1);
f=@(x,y,z)0*x+1;
g=@(x,y,z) 0*x;

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

%% Closest point representations of the solid ball
 
%closest points of the grid points
cpf = @(x,y,z) cpBall(x, y, z, 1, [0 0 0]);
[cpx, cpy, cpz, dist, bdy] = cpf(xx, yy, zz);
  
%find the band
band = find(abs(dist) <= bw*dx);

% store closest points in the band;
cpx = cpx(band); cpy = cpy(band); cpz = cpz(band);
x = xx(band); y = yy(band); z = zz(band);
bdy = bdy(band);
  
% Compute the modified closest points 
[cpbarx, cpbary, cpbarz] = cpbar_3d(x, y, z, cpf);


  
% 3D scatter plot of the cp points
figure(1);clf;
scatter3(cpx, cpy, cpz, 40, 'filled');
hold on;
axis equal
xlabel('x'); ylabel('y'); zlabel('z');
title ('Closest point representation of the ball');
 
%% Build the interpolating matrices and the elliptic operator
% Build the interpolating matrices
E1 = interp3_matrix(x1d, y1d, z1d, cpx, cpy, cpz, 1, band);
E = interp3_matrix(x1d, y1d, z1d, cpbarx, cpbary, cpbarz, p, band);

E(bdy,:) = -E(bdy,:);

% Create Laplacian matrix
L = laplacian_3d_matrix(x1d,y1d,z1d,order, band);
 
% Penalty parameter
gamma = 2*dim/(dx^2) ;
I = speye(size(L)) ;
 
% Modified elliptic operator
M = E1*L - gamma*(I-E);
 
%% Build RHS
rhs = f(cpx, cpy, cpz);
    
% Trilinear interpolation from the grid points to the closest points of the
% ghost_points 
E_bdy = interp3_matrix(x1d, y1d, z1d, cpx(bdy), cpy(bdy), cpz(bdy), 1, band);

% Impose the boundary condition
% rhs(bdy) = E_bdy*g(x,y,z);
bc_dir=zeros(size(x));
bc_dir(bdy)= E_bdy*g(x,y,z);

rhs=rhs - gamma*(2*bc_dir);

% Solve using backslash
tic;
u = M\rhs;
toc

%% Plot of the solution
% Restrict the solution within the ball
 u_in = u(~ (bdy));

% Plot the solution
figure(2);clf;
scatter3(x(~bdy),y(~bdy),z(~bdy),40,u_in, 'filled');
axis equal
xlabel('x'); ylabel('y'); zlabel('z');
title ('Solution of the poisson equation in the ball');
  
%% Error analysis
% The exact solution within the ball
ue_in = uexact(x(~(bdy)), y(~(bdy)), z(~(bdy))); 

% Compute the error in inf norm
err = max(abs(ue_in - u_in))
relerr=max(abs(ue_in - u_in))/max(abs(ue_in))






% uexact=(1/6)* (x.^2+y.^2+z.^2-1);
% f=1;
% g=0;

% dx         Error       Relative error    Elapsed time
% ------------------------------------------------------
% 0.4        0.0031         0.0208           0.031134s
% 0.2        0.0015         0.0087           0.064637s
% 0.1        3.4451e-04     0.0021           1.000462s
% 0.05       9.7140e-05     5.8284e-04       35.726918s






  



