%% Solving an anisotropic elliptic equation on a composite domain in 3D

% The PDE is
%  $$ div [A(y)grad (u(y))] = f$$ in \Omega  
%  $$u=g$$       on\partial\Omega 
%  where \Omega is the composite 3D domain

% The method of manufactured solutions is used to construct
% an exact solution and thus determine f and g.

% Include the cp_matrices folder (edit as appropriate)
addpath('../cp_matrices');
% Add functions for finding the closest points
addpath('../surfaces');
%% Build the mesh grid 
dx = 0.1;

% make vectors of x, y, positions of the grid
x1d = ((-2-6*dx):dx:(2+6*dx))';
y1d = ((-1.5-6*dx):dx:(2+6*dx))';
z1d = ((-2-6*dx):dx:(2+6*dx))';
[xx, yy, zz] = meshgrid(x1d, y1d, z1d);
%% Manufactured solution
uexact= @(x,y,z)x.^2+y.^2+z.^2;
f=@(x,y,z) 10*x.^2+6;
g=uexact;   %A=[x.^2+1, 0 , 0; 0, x.^2+1, 0; 0, 0, x.^2+1]

%% Choose the input 3D objects
which = 6;
% closest point functions
cpf = {};
switch which
  case 1
    cpf{1} = @(x,y,z) cpBall(x, y, z, 0.8, [-0.2 0 0]);
    cpf{2} = @(x,y,z) cpSolidCylinder(x, y, z, [-1, 1], 0.5, [0.5 0]);
  case 2
    cpf{1} = @(x,y,z) cpBall(x, y, z, 0.5, [0.5 0 0.4]);
    cpf{2} = @(x,y,z) cpBall(x, y, z, 0.6, [-0.4 0 -0.3]);
    cpf{3} = @(x,y,z) cpSolidCylinder(x, y, z, [-1, 1], 0.5);
    case 3
      cpf{1} = @(x,y,z) cpBall(x, y, z, 0.5, [0.5 0 0.8]);
      cpf{2} = @(x,y,z) cpBall(x, y, z, 0.6, [-0.4 0 -0.9]);
      cpf{3} = @(x,y,z) cpSolidCylinder(x, y, z, [-1, 0.1], 0.5, [0 0]);
      cpf{4} = @(x,y,z) cpSolidCylinder(x, y, z, [-0.1, 1], 0.5, [0.2 0]);
    case 4
      cpf{1} = @(x,y,z) cpBall(x, y, z, 0.5, [0.5 0 0.4]);
      cpf{2} = @(x,y,z) cpBall(x, y, z, 0.6, [-0.4 0 -0.3]);
      cpf{3} = @(x,y,z) cpSolidCylinder(x, y, z, [-1, 0.1], 0.5);
      cpf{4} = @(x,y,z) cpSolidCylinder(x, y, z, [-0.1, 1], 0.5, [0 0]);
    case 5
      cpf{1} = @(x,y,z) cpBall(x, y, z, 0.5, [-0.5 0 0.4]);
      cpf{2} = @(x,y,z) cpBall(x, y, z, 0.5, [0  0.5 0.4]);
      
    case 6% Mickey Mouse's Face
    cpf{1}=@(x,y,z)cpBall(x,y,z,0.8,[-1,0.4,1]);  
    cpf{2}=@(x,y,z)cpBall(x,y,z,1,[0,0.4,0]);
    cpf{3}=@(x,y,z)cpBall(x,y,z,0.8,[1,0.4,1]);
end

%% Banding
dim = 3;
order = 2;
p = 3;
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));
%% Loop over each object and build banded grids and operators
N = length(cpf);
figure(1);clf;
for i=1:N
  [cpx, cpy, cpz, dist,bdy] = cpf{i}(xx, yy, zz);
  
   band = find(abs(dist) <= bw*dx);
  % store closest points in the band;
  cpx = cpx(band); cpy = cpy(band); cpz = cpz(band);
  x = xx(band); y = yy(band); z = zz(band);
  [cpbarx, cpbary, cpbarz] = cpbar_3d(x, y, z, cpf{i});
  bdy = bdy(band);
  
%   scatter3(cpx(~bdy), cpy(~bdy), cpz(~bdy), 40, 'k','filled');
scatter3(cpx, cpy, cpz, 40, 'k','filled');
  hold on;
  az = -20;
  el = 10;
  view(az, el);

  E1 = interp3_matrix(x1d, y1d, z1d, cpx, cpy, cpz, 1, band);
  E = interp3_matrix(x1d, y1d, z1d, cpx, cpy, cpz, p, band);

%%  Build the differential operator and Interpolation matrices
 rho=@(x,y,z) x.^2+1;
 P=diag(sparse(rho(x,y,z)));
 
%  Build the 3D averaging matrices
[Axb,Axf,Ayb,Ayf,Azb,Azf] = average_upw1_3d_matrices(x1d,y1d,z1d,band);

%  Build the 3D difference matrices
[Dxb,Dxf,Dyb,Dyf,Dzb,Dzf]=firstderiv_upw1_3d_matrices(x1d,y1d,z1d,band);

% Build the differential operator
L=Dxb*(Axf*P*Dxf)+Dyb*(Ayf*P*Dyf)+Dzb*(Azf*P*Dzf);

% Penalty parameter
gamma = 2*dim/(dx^2) ;
I = speye(size(L)) ;
 
% Modified elliptic operator
M = E1*L - gamma*(I-E);

 %% Deal with the boundary
% Compute the modified closest points which will be used to deal with boundary
[cpbarx, cpbary,cpbarz] = cpbar_3d(x, y,z, cpf{i});

% Build the interpolating matrix for ghost points
E_bdy=interp3_matrix(x1d,y1d,z1d,cpbarx(bdy), cpbary(bdy),cpbarz(bdy), p ,band);

% Modidy the entries of the operator matrix
M(bdy,:)=(I(bdy,:)+E_bdy)/2;

%% Build RHS
rhs = f(cpx, cpy,cpz);
    
% Trilinear interpolation from the grid points to the closest points of the
% ghost_points 
E_ghost = interp3_matrix(x1d, y1d,z1d, cpx(bdy), cpy(bdy),cpz(bdy), 1, band);

% Impose the boundary condition
rhs(bdy) = E_ghost*g(x,y,z); 

d{i}.bdy = bdy;
d{i}.cpx = cpx;
d{i}.cpy = cpy;
d{i}.cpz = cpz;
d{i}.cpbarx=cpbarx;
d{i}.cpbary=cpbary;
d{i}.cpbarz=cpbarz;
d{i}.x = x;
d{i}.y = y;
d{i}.z = z;
d{i}.M = M;
d{i}.I = I;
d{i}.rhs = rhs;
d{i}.band=band;
end

%% Initial guess u0
uold = {};
for i=1:N
    s= d{i}; 
    uold{i} = cos(10*s.x) + sin(12*s.y)+s.z;
end

% 3D plot of the initial guess
figure(2); clf;
for i=1:N
    s= d{i}; 
    scatter3(s.cpx, s.cpy,s.cpz, 40, uold{i}, 'filled');
    hold on;
    title('Plot of u^{0}')
    xlabel('x'); ylabel('y');zlabel('z'); 
    az = -20;
    el = 10;
    view(az, el);
    drawnow   
end


%% outer iteration (swartz?)
for J = 1:20
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
    [C, Io, Is] = intersect(o.band(~o.bdy), s.band(s.bdy));
    temp = find(~o.bdy);  Io = temp(Io);
    temp = find(s.bdy);  Is = temp(Is);
    assert (all (C == o.band(Io)));
    assert (all (C == s.band(Is)));
 
    E_bdy_in=interp3_matrix(x1d, y1d,z1d, s.cpbarx(Is), s.cpbary(Is),s.cpbarz(Is), 1, s.band);
    s.M(Is,:)=(s.I(Is,:)+E_bdy_in)/2;
% Trilinear interpolation from the grid points to the closest points of the
% ghost_points (with respect to the other sub-domains)
 E_ghost_in= interp3_matrix(x1d, y1d,z1d, s.cpx(Is), s.cpy(Is), s.cpz(Is), 1,o.band);
 s.rhs(Is) = E_ghost_in*uold{j};
    
    % Plot of the ghost points
    % figure(10+i); clf; 
    % scatter3(s.cpx(Is),s.cpy(Is),s.cpz(Is), 40, 'filled')
      
    %figure(10+i); clf; I = s.bdy
    %scatter3(xx(I),yy(I),zz(I), 40, zz(I)+20, 'filled')
    %axis equal; hold on; pause
  end

  tic
  u{i}= (s.M)\s.rhs;
  toc
end

err = [J];
for i=1:N
  s = d{i};
  ue_in = uexact(s.x(~(s.bdy)),s.y(~(s.bdy)),s.z(~(s.bdy)));
  u_in = u{i}(~(s.bdy));
  err(i+1) = max(abs(ue_in - u_in));
  relerr=max(abs(ue_in-u_in))/(max(abs(ue_in)))
end
err
uold = u;

figure(3); clf;
for i=1:N
  s = d{i};
  scatter3(s.cpx,s.cpy, s.cpz, 40, uold{i}, 'filled')
  hold on;
  axis equal
  title(['3D plot of u^{' num2str(J) '} for CP Additive Schwarz Method'])
  xlabel('x'); ylabel('y');zlabel('z');
  drawnow
  az = -20;
  el = 10;
  view(az, el);
end
end