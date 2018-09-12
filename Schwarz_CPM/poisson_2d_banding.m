function u=poisson_2d_banding(f,g,dx,x1d,y1d,cpf,varargin)
% poisson_2d_banding solves 2d Poisson equation using the Closest Point
% Method with banding
%  u=poisson_2d_banding(f,g, x1d,y1d,cpf,varargin) solves the two dimensional equation
%   -\Delta u=f on the domain Omega with the cloeset point representations [cpx cpy] 
%    with Dirichlet boundary conditions u=g 
%% The Closest Point Representations
[x,y]=meshgrid(x1d,y1d);
% make into vectors
x=x(:); y=y(:);
[cpx, cpy, dist,bdy]=cpf(x,y,varargin{:});
%% Banding: do calculation in a narrow band around the circle
dim = 2;
p = 3;
order=2;
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
band = find(abs(dist) < bw*dx);

% store closest pionts in the band;
cpx = cpx(band); cpy = cpy(band); dist= dist(band);; 
x = x(band); y = y(band);
bdy = bdy(band);
% Find the modified closest point representation of each grid point
% "cpbar" [Macdonald, Brandman, Ruuth 2011]:cpbar(x):=cp(2*cp(x)-x)
[cpbar_x, cpbar_y]=cpf(2*cpx-x, ...
2*cpy-y,varargin{:});
cpbar_x=cpbar_x(:);
cpbar_y=cpbar_y(:);

%% Build the Laplacian operator
% For the second-order Laplace-Beltrami operator, we need p>=q+1 to achieve
% the full order of accuracy of the underlying finite difference scheme
% (say q) [Macdonald & Ruuth 2009]
disp('building laplacian and interp matrices');
L = laplacian_2d_matrix(x1d,y1d, order, band, band);
E = interp2_matrix(x1d,y1d,cpx, cpy, p,band);
E1 = interp2_matrix(x1d,y1d,cpx,cpy,1,band);
I = speye(size(L));
% Penalty parameter 
gamma=2*dim/dx^2;

% Modified Laplacian operator [Macdonald, von Glehn and März,2014]
L= E1*I*L-gamma*(I - E); 
M=-L;

%% Building matrices to deal with dirichlet boundary conditions 
disp('buidling matrices to deal with boundary conditions ... ');
E_bar=interp2_matrix(x1d,y1d,cpbar_x(bdy), cpbar_y(bdy), p,band);
M_bdy = (I(bdy,:) + E_bar)/2;
M(bdy,:) = M_bdy;
disp('done');
%  Bilinear interpolation from the grid points to the closest points 
E_bdy=interp2_matrix(x1d,y1d,cpx(bdy),cpy(bdy),1,band);
%% Solve the Poisson Equation
f=f(cpx,cpy);
f(bdy)= E_bdy*g(cpx,cpy);
u=M\f;
%% Building the Plotting matrices
disp('building plotting matrices');

Eplot= interp2_matrix(x1d,y1d,x(~bdy),y(~bdy),p,band);
% Plot the solution to the whole domain
figure(1);
plot2d_compdomain2(Eplot*u,x(~bdy),y(~bdy),dx,dx,1);
end
