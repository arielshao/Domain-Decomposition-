%% Solving an anisotropic elliptic equation on a union of 2D-objects using Schwarz iterations 

%  This example solves an anisotropic elliptic equation on a union of
%  2D-objects
%  with dirichlet boundary condition.

%  The PDE is
%  $$ div [A(y)grad (u(y))] = f$$ in \Omega  
%  $$u=g$$       on\partial\Omega 
%  where \Omega is the composite 2D domain

% The method of manufactured solutions is used to construct
% an exact solution and thus determine f and g.


% Include the cp_matrices folder (edit as appropriate)
addpath('../cp_matrices');
% add functions for finding the closest points
addpath('../surfaces');

%% Manufactured solution
uexact= @(x,y)x.^2+y.^2;
f=@(x,y) 8*x.^2+4;
g=uexact;   %A=[x.^2+1, 0; 0, x.^2+1]

% uexact=@(x,y)x.^2+y.^2;
% f=@(x,y) 6*(x.^2+y.^2)+4+2*y;
% g=uexact;  % A=[x.^2+1, x; x, y.^2+1]
%% Build the mesh grid 
dx=0.05; % grid size
x1d=(-2-5*dx:dx:2+5*dx)';
y1d=(-2-5*dx:dx:2+5*dx)';
% Generate meshgrid
[xx,yy]=meshgrid(x1d,y1d);
% % Make into vectors
% X=xx(:); Y=yy(:);
% % Plot the grid points
% plot(X, Y, '.', 'color', 0.75*[1 1 1]);
% hold on;
%% Choose the input 2D objects
R=1;
r=0.8; 
which =5;
% closest point functions
cpf = {};
switch which
  case 1
    cpf{1}=@(x,y)cpAnnulus2(x,y,r,R,[-3,0]);  
    cpf{2}=@(x,y)cpAnnulus2(x,y,r,R,[-1.5,-0.5]);
    cpf{3}=@(x,y)cpAnnulus2(x,y,r,R,[0,0]);
    cpf{4}=@(x,y)cpAnnulus2(x,y,r,R,[1.5,-0.5]);
    cpf{5}=@(x,y)cpAnnulus2(x,y,r,R,[3,0]);
  case 2
    cpf{1}=@(x,y)cpAnnulus2(x,y,r,R,[-3,0]);  
    cpf{2}=@(x,y)cpAnnulus2(x,y,r,R,[-1.5,-0.5]);
    cpf{3}=@(x,y)cpAnnulus2(x,y,r,R,[0,0]);
    cpf{4}=@(x,y)cpAnnulus2(x,y,r,R,[1.5,-0.5]);
  case 3
    cpf{1}=@(x,y)cpAnnulus2(x,y,r,R,[-3,0]);  
    cpf{2}=@(x,y)cpAnnulus2(x,y,r,R,[-1.5,-0.5]);
    cpf{3}=@(x,y)cpAnnulus2(x,y,r,R,[0,0]);
  case 4
     cpf{1}=@(x,y)cpAnnulus2(x,y,r,R,[-3,0]);  
     cpf{2}=@(x,y)cpAnnulus2(x,y,r,R,[-1.5,-0.5]);       
  case 5 % Mickey Mouse's Face
    cpf{1}=@(x,y)cpDisc(x,y,0.8,[-1,1]);  
    cpf{2}=@(x,y)cpDisc(x,y,1,[0,0]);
    cpf{3}=@(x,y)cpDisc(x,y,0.8,[1,1]);
end

%% Banding
dim = 2;
order = 2;
p = 3;
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));

%% Loop over each object and build banded grids and operators
N = length(cpf);
figure(1);clf;
for i=1:N
  [cpx, cpy, dist, bdy] = cpf{i}(xx, yy);
   band = find(abs(dist) <= bw*dx);

  % store closest points in the band;
  cpx = cpx(band); cpy = cpy(band);
  x = xx(band); y = yy(band);
  bdy = bdy(band);
  
  
 
scatter(cpx(~bdy), cpy(~bdy), 40, 'k','filled'); hold on;
xlabel('x'); ylabel('y');
 xlim([-2, 2]);
 ylim([-1.5, 2]);
  E1 = interp2_matrix(x1d, y1d,cpx, cpy,1,band);
  E = interp2_matrix(x1d, y1d,cpx, cpy,p, band);


%  Build the differential operator and Interpolation matrices
%  Generate the function-valued matrix A  
  A={@(x,y)x.^2+1, @(x,y) 0; @(x,y) 0,@(x,y) x.^2+1}
  
%  A={@(x,y)x.^2+1, @(x,y) x; @(x,y) x,@(x,y) y.^2+1}

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

% P11=diag(sparse(a11(cpx,cpy)));
% P12=diag(sparse(a12(cpx,cpy)));
% P21=diag(sparse(a21(cpx,cpy)));
% P22=diag(sparse(a22(cpx,cpy)));

%  Build the 2D averaging matrices
[Axb,Axf,Ayb,Ayf]=average_upw1_2d_matrices(x1d,y1d,band);
% Axf(1,1)=1; Axb(1,1)=1;Ayf(1,1)=1; Ayb(1,1)=1;
% Axf(end,end)=1; Axb(end,end)=1; Ayf(end,end)=1; Ayb(end,end)=1; 

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
% Compute the modified closest points which will be used to deal with boundary
[cpbarx, cpbary] = cpbar_2d(x, y, cpf{i});

% Build the interpolating matrix for ghost points
E_bdy=interp2_matrix(x1d,y1d,cpbarx(bdy), cpbary(bdy), p ,band);

% Modidy the entries of the operator matrix
M(bdy,:)=(I(bdy,:)+E_bdy)/2;

%% Build RHS
rhs = f(cpx, cpy);
    
% bilinear interpolation from the grid points to the closest points of the
% ghost_points 
E_ghost = interp2_matrix(x1d, y1d,cpx(bdy), cpy(bdy), 1, band);


% Impose the boundary condition
rhs(bdy) = E_ghost*g(x,y); 

d{i}.bdy = bdy;
d{i}.cpx = cpx;
d{i}.cpy = cpy;
d{i}.cpbarx=cpbarx;
d{i}.cpbary=cpbary;
d{i}.x = x;
d{i}.y = y;
d{i}.M = M;
d{i}.I = I;
d{i}.rhs = rhs;
d{i}.band=band;
end

%% Initial guess u0
uold = {};
for i=1:N
    s= d{i}; 
    uold{i} = cos(10*s.x) + sin(12*s.y);
end

% 2D plot of the initial guess
figure(N+1); clf;

for i=1:N
    s= d{i}; 
    scatter(s.cpx(~s.bdy), s.cpy(~s.bdy), 40, uold{i}(~s.bdy), 'filled');
    xlim([-2, 2]);
    ylim([-1.5, 2]);
    hold on;
    axis equal
    title('Plot of u^{0}')
    xlabel('x'); ylabel('y'); 
    drawnow
end


%% outer iteration (swartz?)

Nsteps=30;
for J = 1:Nsteps;
u = {};
for i=1:N
  s = d{i};

  %% build a randomly ordered list of the other shapes
  jset = setdiff(randperm(N), i);
  for j=jset
      
      o = d{j};
    %% now "s" is the current object we're solving on and "o" is the
    % other object.  Locate ghost points for s that are in the
    % interior of o.
 
     % find the ghost points
    
    [C, Io, Is] = intersect(o.band(~o.bdy), s.band(s.bdy));
    temp = find(~o.bdy);  Io = temp(Io);
    temp = find(s.bdy);  Is = temp(Is);
    assert (all (C == o.band(Io)));
    assert (all (C == s.band(Is)));
  
    E_bdy_in=interp2_matrix(x1d, y1d, s.cpbarx(Is), s.cpbary(Is), 1, s.band);
    s.M(Is,:)=(s.I(Is,:)+E_bdy_in)/2;
% Trilinear interpolation from the grid points to the closest points of the
% ghost_points (with respect to the other sub-domains)

%    E_ghost_in= interp2_matrix(x1d, y1d,  s.cpx(Is), s.cpy(Is), 1, o.band);
 E_ghost_in= interp2_matrix(x1d, y1d, s.cpx(Is), s.cpy(Is), 1,o.band);
 s.rhs(Is) = E_ghost_in*uold{j};
 end
 tic
 u{i}= (s.M)\s.rhs;
 toc
end
err = [J];
for i=1:N
  s = d{i};
  ue_in = uexact(s.x(~(s.bdy)),s.y(~(s.bdy)));
  u_in = u{i}(~(s.bdy));
  err(i+1) = max(abs(ue_in - u_in));
  relerr=max(abs(ue_in-u_in))/(max(abs(ue_in)))
end
err
uold = u;
 for i=1:N
  s = d{i};
  scatter(s.cpx,s.cpy,40, uold{i}, 'filled')
  hold on;
  xlim([-2, 2]);
  ylim([-1.5, 2]);
  title(['2D plot of u^{' num2str(J) '} for CP Additive Schwarz Method'])
  xlabel('x'); ylabel('y');
  drawnow
 end
end