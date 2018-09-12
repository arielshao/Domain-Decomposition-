%% Solving an anisotropic elliptic equation on an elliptical-shape disc with Dirichlet B.C

%  This example solves an anisotropic elliptic equation on an elliptical-shaped disc 
%  with dirichlet boundary condition.

%  The PDE is
%  $$ div [A(y)grad (u(y))] = f$$ in \Omega  
%  $$u=g$$       on\partial\Omega 
%  where \Omega is the elliptical-shape disc


% This code deals with the boundary condition using
% linear extrapolation, that is,
% u(x_g)=2u(cp(x_g))-u(cpbar(x_g)) where cpbar=cp(2cp(x)-x)


% The method of manufactured solutions is used to construct
% an exact solution and thus determine f and g.

% adjust as appropriate
addpath('../cp_matrices');
addpath('../surfaces');

%% Manufactured solution
uexact= @(x,y)x.^2+y.^2;
f=@(x,y) 8*x.^2+4;
g=uexact;   %A=[x.^2+1, 0; 0, x.^2+1]


% uexact=@(x,y)x.^2+y.^2;
% f=@(x,y) 6*(x.^2+y.^2)+4+2*y;
% g=uexact;  % A=[x.^2+1, x; x, y.^2+1]
% 


%% Build the mesh grid 

dx = 0.005;  %grid size

% Construct the grid in the embedding space
x1d = ((-1-6*dx):dx:(1+6*dx))';
y1d = ((-1-6*dx):dx:(1+6*dx))';

[xx, yy] = meshgrid(x1d, y1d);

%% Banding (No banding in this example)
dim = 2;
order = 2;
p = 3;
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));
%% Plot the ellipse
LW='LineWidth';
a=1; % Horizontal radius
b=0.5; % Vertical radius 
centrex=0;  % [x0,y0] ellipse centre coordinates
centrey=0;
t=-pi:0.01:pi;
xcoor=centrex+a*cos(t);
ycoor=centrey+b*sin(t);
figure(1);clf
plot(xcoor,ycoor,'k-');
xlim([-1.2 1.2]);
ylim([-1.2 1.2]);
hold on;

%% Closest point representations of the elliptical disc
 %closest points of the grid points
cpf = @(x,y) cpEllipticalDisc(x, y, a, b, [0 0]);
[cpx, cpy, dist, bdy] = cpf(xx, yy);
  
%find the band
band = find(abs(dist) <=bw*dx);

%store closest points in the band;
cpx = cpx(band); cpy = cpy(band); 
x = xx(band); y = yy(band); 
bdy = bdy(band);

% Compute the modified closest points which will be used to deal with boundary
[cpbarx, cpbary] = cpbar_2d(x,y,cpf);
   
% Plot the grid points
plot(x, y, '.', 'color', 0.75*[1 1 1]); 
hold on;

% 2D scatter plot of the cp points
scatter(cpx, cpy, 40, 'filled');
axis equal
xlabel('x'); ylabel('y'); 
title ('Closest point representation of the elliptical disc');

%% Build the interpolating matrices and the elliptic operator
% Build the interpolating matrices
E1 = interp2_matrix(x1d, y1d, cpx, cpy, 1, band);
E = interp2_matrix(x1d, y1d,cpbarx, cpbary, p, band);
% Build the differential operator and Interpolation matrices
% Generate the function-valued matrix A  
   A={@(x,y)x.^2+1, @(x,y) 0; @(x,y) 0,@(x,y) x.^2+1}
% A={@(x,y)x.^2+1, @(x,y) x; @(x,y) x,@(x,y) y.^2+1}

 %We need our matrix A=[a11(x,y) a12(x,y) ; a21(x,y) a22(x,y)] to be positive-definite
 % Note that if we take a11=a22=x^2, a12=a21=0,then the matrix is zero at y-axis 
 % and this will cause discontinuity of the plot along y-axis
 
 % Change the cells to functions
 a11=cell2mat(A(1,1));
 a12=cell2mat(A(1,2));
 a21=cell2mat(A(2,1));
 a22=cell2mat(A(2,2));
 
% Build the diagonal matrices using the functions from matrix A
P11=diag(sparse(a11(x,y)));
P12=diag(sparse(a12(x,y)));
P21=diag(sparse(a21(x,y)));
P22=diag(sparse(a22(x,y)));

%  Build the 2D averaging matrices
[Axb,Axf,Ayb,Ayf]=average_upw1_2d_matrices(x1d,y1d,band);
Axf(1,1)=1; Axb(1,1)=1;Ayf(1,1)=1; Ayb(1,1)=1;
Axf(end,end)=1; Axb(end,end)=1; Ayf(end,end)=1; Ayb(end,end)=1; 

[Dxb,Dxf,Dyb,Dyf]=firstderiv_upw1_2d_matrices(x1d,y1d,band);
% Build the differential operator
% L=Dxb*(Axf*P11*Dxf)+Dyb*(Ayf*P22*Dyf);
 L=Dxb*(Axf*P11*Dxf)+Dyb*(Ayf*P22*Dyf)+Dyb*(Axf*P21*Dxf)+Dxb*(Ayf*P12*Dyf);

 
% Penalty parameter
gamma = 2*dim/(dx^2) ;
I = speye(size(L)) ;
 
% Modified elliptic operator
M = E1*L - gamma*(I-E);

%% Deal with the boundary
% % Compute the modified closest points which will be used to deal with boundary
[cpbarx, cpbary] = cpbar_2d(x,y,cpf);

% Build the interpolating matrix for ghost points
E_bdy=interp2_matrix(x1d,y1d,cpbarx(bdy), cpbary(bdy),p, band);

% Modidy the entries of the operator matrix
M(bdy,:)=(I(bdy,:)+E_bdy)/2;
%% Build RHS
rhs = f(cpx, cpy);
    
% Bilinear interpolation from the grid points to the closest points of the
% ghost_points 
E_ghost = interp2_matrix(x1d, y1d, cpx(bdy), cpy(bdy), 1,band);
% Impose the boundary condition
rhs(bdy) = E_ghost*g(x,y);    

% Solve using backslash
tic;
u = M\rhs;
toc

%% Plot of the solution
% Restrict the solution within the elliptical-shape disc
 u_in = u(~(bdy));

% Plot the solution
figure(2);clf;
% scatter(xx(~bdy),yy(~bdy), 40,u_in,'filled');
plot2d_compdomain2(u_in,x(~bdy),y(~bdy),dx,dx,2);
hold on;
plot(xcoor,ycoor,'k-');
xlim([-1.2 1.2]);
ylim([-1.2 1.2]);
xlabel('x'); ylabel('y'); 
title ('Solution of the anisotropic elliptic equation on the elliptical-shape disc');
%   
%% Error analysis
% The exact solution within the elliptical-shape disc
ue_in = uexact(x(~(bdy)), y(~(bdy))); 

% Compute the error in inf norm
err=max(abs(ue_in-u_in))
relerr=max(abs(ue_in-u_in))/max(abs(ue_in))

%% Observations and Analysis:

% If %A=[x.^2+1, 0; 0, x.^2+1]
% uexact= @(x,y)x.^2+y.^2; f=@(x,y) 8*x.^2+4; g=uexact;   

% dx         Error       Relative error    Elapsed time
% ------------------------------------------------------
% 0.02       0.0028        0.0029          0.035635s
% 0.01       8.3200e-04    8.4467e-04      0.221556s
% 0.005      2.6243e-04    2.6454e-04      2.892799s

% Though the error is small with small grid size, it does not seem to
% converge

%  If A=[x.^2+1, x; x, y.^2+1]
% uexact=@(x,y)x.^2+y.^2; f=@(x,y) 6*(x.^2+y.^2)+4+2*y; g=uexact;  
% The solution blows up






