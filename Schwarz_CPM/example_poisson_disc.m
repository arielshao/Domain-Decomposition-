%% Solving a Poisson's problem on a disc with Dirichlet B.C
%  This example solves the poisson's equation on a disc 
%  with dirichlet boundary condition.

%  The PDE is
%  $$lap u = f$$ in \Omega  
%  $$u=g$$       on\partial\Omega 
%  where \Omega is the disc


% This code deals with the boundary condition using
% linear extrapolation, that is,
% u(x_g)=2u(cp(x_g))-u(cpbar(x_g)) where cpbar=cp(2cp(x)-x)


% The method of manufactured solutions is used to construct
% an exact solution and thus determine f and g.

% adjust as appropriate
addpath('../cp_matrices');
addpath('../surfaces');

 %% Manufactured solution
uexact=@(x,y) (1/4)* (x.^2+y.^2-1);
f=@(x,y,z)0*x+1;
g=@(x,y,z) 0*x;

%% Build the mesh grid 

dx = 0.05;  %grid size

% Construct the grid in the embedding space
x1d = ((-1-6*dx):dx:(1+6*dx))';
y1d = ((-1-6*dx):dx:(1+6*dx))';

[xx, yy] = meshgrid(x1d, y1d);


%% Banding
dim = 2;
order = 2;
p = 3;
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));
%% Closest point representations of the disc
 
%closest points of the grid points
cpf = @(x,y,z) cpDisc(x, y, 1, [0 0]);
[cpx, cpy, dist, bdy] = cpf(xx, yy);
  
%find the band
band = find(abs(dist) <= bw*dx);

% store closest points in the band;
cpx = cpx(band); cpy = cpy(band); 
x = xx(band); y = yy(band); 
bdy = bdy(band);
   
figure(1);clf;
% Plot the grid points
plot(x, y, '.', 'color', 0.75*[1 1 1]); 
hold on;

% 2D scatter plot of the cp points
scatter(cpx, cpy, 40, 'filled');
% Plot the circle
thetas=linspace(0,2*pi,100)';
r=ones(size(thetas));
[xcir,ycir]=pol2cart(thetas,r);
plot(xcir(:), ycir(:), 'k-', 'LineWidth', 2);
hold on;
xlim([-1.2, 1.2]);
ylim([-1.2, 1.2]);
axis equal
xlabel('x'); ylabel('y'); 
title ('Closest point representation of the disc');

%% Build the interpolating matrices and the elliptic operator
% Build the interpolating matrices
E1 = interp2_matrix(x1d, y1d, cpx, cpy, 1, band);
E = interp2_matrix(x1d, y1d,cpx, cpy, p, band);
  
% Create Laplacian matrix
L = laplacian_2d_matrix(x1d,y1d,order, band);
 
% Penalty parameter
gamma = 2*dim/(dx^2) ;
I = speye(size(L)) ;
 
% Modified elliptic operator
M = E1*L - gamma*(I-E);


%% Deal with the boundary
% Compute the modified closest points which will be used to deal with boundary
[cpbarx, cpbary] = cpbar_2d(x, y, cpf);

% Build the interpolating matrix for ghost points
E_bdy=interp2_matrix(x1d,y1d,cpbarx(bdy), cpbary(bdy), p,band);

% Modidy the entries of the operator matrix
M(bdy,:)=(I(bdy,:)+E_bdy)/2;
%% Build RHS
rhs = f(cpx, cpy);
    
% Bilinear interpolation from the grid points to the closest points of the
% ghost_points 
E_ghost = interp2_matrix(x1d, y1d, cpx(bdy), cpy(bdy),  1, band);

% Impose the boundary condition
rhs(bdy) = E_ghost*g(x,y);    

% Solve using backslash
tic;
u = M\rhs;
toc

%% Plot of the solution
% Restrict the solution within the disc
 u_in = u(~ (bdy));

% Plot the solution
figure(2);clf;
scatter(x(~bdy),y(~bdy), 40,u_in,'filled');

% plot2d_compdomain2(u_in,x(~bdy),y(~bdy),dx,dx,2);
 hold on;
% Plot the circle
plot(xcir(:), ycir(:), 'k-', 'LineWidth', 2);
xlim([-1.2, 1.2]);
ylim([-1.2, 1.2]);
axis equal
xlabel('x'); ylabel('y'); 
title ('Solution of the poisson equation on the disc');
  
%% Error analysis
% The exact solution within the disc
ue_in = uexact(x(~(bdy)), y(~(bdy))); 

% Compute the error in inf norm
err = max(abs(ue_in - u_in))
relerr=max(abs(ue_in - u_in))/max(abs(ue_in))




