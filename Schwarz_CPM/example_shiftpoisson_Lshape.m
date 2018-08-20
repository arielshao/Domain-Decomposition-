%% Solving a shift Poisson's problem on a Lshape domain with Neumann B.C
%  This example solves the shift poisson's equation on a Lshape domain
%  with neumann boundary condition.

%  The PDE is
%  $$-\lap u+u = f$$ in \Omega  
%  $$ \frac{\partial u}{\partial n}=g$$     on\partial\Omega 
%  where \Omega is the Lshape domain

% Note that cp extension gives zero neumann b.c.

% The method of manufactured solutions is used to construct
% an exact solution and thus determine f and g.

% adjust as appropriate
addpath('../cp_matrices');
addpath('../surfaces');

%% Manufactured solution
uexact=@(x,y) 0*x+1;
f=@(x,y,z)0*x+1;
% g=@(x,y,z) 0*x;  
%% Build the mesh grid 
dx = 0.05;  %grid size

% Construct the grid in the embedding space
x1d = ((-1.5-6*dx):dx:(1.5+6*dx))';
y1d = ((-1.5-6*dx):dx:(1.5+6*dx))';
[xx, yy] = meshgrid(x1d, y1d);

%% Banding
dim = 2;
order = 2;
p = 3;
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));
%% Closest point representations of the Lshape domain
%closest points of the grid points
cpf = @(x,y,z) cpLshapeDomain(x, y, 2,1, [0 0]);
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
hold on;
% Plot the L-shape domain
% Plot the lower-rectangular domain 
xlow1=-1;xlow2=1;
ylow1=-1;ylow2=0;
xlow = [xlow1, xlow2, xlow2, xlow1, xlow1];
ylow = [ylow1, ylow1, ylow2, ylow2, ylow1];
plot(xlow, ylow, 'b-', 'LineWidth', 2);
xlim([-1.5 1.5]);
ylim([-1.5 1.5]);
hold on;
% Plot the upper-rectangular domain 
xup1=0;xup2=1;
yup1=-1;yup2=1;
xup = [xup1, xup2, xup2, xup1, xup1];
yup = [yup1, yup1, yup2, yup2, yup1];
plot(xup, yup, 'b-', 'LineWidth', 2);
axis equal
xlabel('x'); ylabel('y'); 
title ('Closest point representation of the Lshape domain');

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
M =-( E1*L - gamma*(I-E))+I;

%% Build RHS
rhs = f(cpx, cpy);   

% Solve using backslash
tic;
u = M\rhs;
toc

%% Plot of the solution
% Restrict the solution within the Lshape domain
 u_in = u(~ (bdy));

% Plot the solution
figure(2);clf;
scatter(x(~bdy),y(~bdy), 40,u_in,'filled');

% plot2d_compdomain2(u_in,x(~bdy),y(~bdy),dx,dx,2);
 hold on;
% Plot of the Lshape domain
plot(xup, yup, 'b-', 'LineWidth', 2);
plot(xlow, ylow, 'b-', 'LineWidth', 2);
xlim([-1.5, 1.5]);
ylim([-1.5, 1.5]);
axis equal
xlabel('x'); ylabel('y'); 
title ('Solution of the poisson equation on the Lshape Domain');
  
%% Error analysis
% The exact solution within the disc
ue_in = uexact(x(~(bdy)), y(~(bdy))); 

% Compute the error in inf norm
err = max(abs(ue_in - u_in))
relerr=max(abs(ue_in - u_in))/max(abs(ue_in))



