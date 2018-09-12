function u=poisson_2d(f,g,dx,x1d,y1d,cpf,varargin)
% poisson_2d solves 2d Poisson equation using the Closest Point Method
%  u=poisson_2d(f,g, x1d,y1d,cpf,varargin) solves the two dimensional equation
%   -\Delta u=f on the domain Omega with the cloeset point representations [cpx cpy] 
%    with Dirichlet boundary conditions u=g 
%% The Closest Point Representations
[x,y]=meshgrid(x1d,y1d);
% make into vectors
x=x(:); y=y(:);
[cpx, cpy, dist1,bdy1]=cpf(x,y,varargin{:});
% Find the modified closest point representation of each grid point
% "cpbar" [Macdonald, Brandman, Ruuth 2011]:cpbar(x):=cp(2*cp(x)-x)
[cpbar_x, cpbar_y,dist2,bdy2]=cpf(2*cpx-x, ...
2*cpy-y,varargin{:});
cpbar_x=cpbar_x(:);
cpbar_y=cpbar_y(:);
%% Build the Laplacian operator
dim = 2;    % dimension
p = 3;      % interpolation degree
% For the second-order Laplace-Beltrami operator, we need p>=q+1 to achieve
% the full order of accuracy of the underlying finite difference scheme
% (say q) [Macdonald & Ruuth 2009]
disp('building laplacian and interp matrices');
[Dxx,Dyy,Dxc,Dyc,Dxb,Dyb,Dxf,Dyf,Dxyc]= diff2d_matrices(x1d,y1d);
L=Dxx+Dyy;
E = interp2_matrix(x1d,y1d,cpx, cpy, p);
E1 = interp2_matrix(x1d,y1d,cpx, cpy,1);
I = speye(size(L));
% Penalty parameter 
gamma=2*dim/dx^2;

% Modified Laplacian operator [Macdonald, von Glehn and März,2014]
L= E1*I*L-gamma*(I - E); 
M=-L;
%% Find the ghost points we need for the Dirichlet boudnary part
ghost_Gamma=(sqrt((cpbar_x-x).^2+(cpbar_y-y).^2)>...
sqrt((x-cpx).^2+(y-cpy).^2));

%% Building matrices to deal with dirichlet boundary conditions 
disp('buidling matrices to deal with boundary conditions ... ');
E_Gamma_bar=interp2_matrix(x1d,y1d,cpbar_x(ghost_Gamma), cpbar_y(ghost_Gamma), p);
M_ghost_Gamma = (I(ghost_Gamma,:) + E_Gamma_bar)/2;
M(ghost_Gamma,:) = M_ghost_Gamma;
disp('done');
%  Bilinear interpolation from the grid points to the closest points 
E_ghost_Gamma=interp2_matrix(x1d,y1d,cpx(ghost_Gamma),cpy(ghost_Gamma),1);
%% Solve the Poisson Equation
f=f(cpx,cpy);
f(ghost_Gamma)= E_ghost_Gamma*g(cpx,cpy);
u=M\f;
%% Building the Plotting matrices
disp('building plotting matrices');
%Get grid points inside the disc
tol=dx/2;
plot_index=(sqrt((cpbar_x-x).^2+(cpbar_y-y).^2)-sqrt((x-cpx).^2+(y-cpy).^2)<tol);
xx=x(plot_index);yy=y(plot_index);
Eplot= interp2_matrix(x1d,y1d,xx,yy,p);
% Plot the solution to the whole domain
figure(1);
plot2d_compdomain2(Eplot*u,xx,yy,dx,dx,1);
end

