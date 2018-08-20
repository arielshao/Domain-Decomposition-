%% Solving a Poisson problem on a disc using Schwarz iterations 

%  This example solves the poisson's equation on a disc
%  with dirichlet boundary condition using Schwarz iterations

%  The PDE is
%  $$lap u = f$$ in \Omega  
%  $$u=g$$       on\partial\Omega 
%  where \Omega is the disc consisting of a union of 2 half discs.

% The method of manufactured solutions is used to construct
% an exact solution and thus determine f and g.


% Include the cp_matrices folder (edit as appropriate)
addpath('../cp_matrices');
% add functions for finding the closest points
addpath('../surfaces');
%% Build the mesh grid 
dx =0.05;

% make vectors of x, y, positions of the grid
x1d = ((-1-6*dx):dx:(1+6*dx))';
y1d = ((-1-6*dx):dx:(1+6*dx))';

[xx, yy] = meshgrid(x1d, y1d);

 %% Manufactured solution
uexact=@(x,y) (1/4)* (x.^2+y.^2-1);
f=@(x,y,z)0*x+1;
g=@(x,y,z) 0*x;

%% Banding
dim = 2;
order = 2;
p = 3;
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));

%% Closest point representations of the 2D subdomains 
cpf = {};
%closest points of the grid points
cpf{1}=@(x,y) cpDisc_left(x,y,1,  [0 0], 0.3);
cpf{2}=@(x,y) cpDisc_right(x,y,1, [0 0], -0.3);

%% Build the elliptic operator and interpolating matrices for each subdomain

N=length(cpf);
for i=1:N

[cpx, cpy, dist, bdy] = cpf{i}(xx, yy);
band = find(abs(dist) <= bw*dx);
 
% store closest points in the band;
cpx= cpx(band); cpy = cpy(band); 
x = xx(band); y = yy(band); 
bdy = bdy(band);

% 2D scatter plot of the cp points
figure(i);clf;
scatter(cpx, cpy, 40, 'filled');
xlim([-1.2, 1.2]);
ylim([-1.2, 1.2]);
hold on;
% Plot the grid points
plot(x, y, '.', 'color', 0.75*[1 1 1]); 
hold on;
axis equal
xlabel('x'); ylabel('y'); 
title (['Closest point representation of the subdomain ' num2str(i)]);

% Build the interpolating matrices
E1 = interp2_matrix(x1d, y1d, cpx, cpy,  1, band);
E = interp2_matrix(x1d, y1d,  cpx, cpy,  p, band);
  
% Create Laplacian matrix
L = laplacian_2d_matrix(x1d,y1d,order, band);
 
% Penalty parameter
gamma = 2*dim/(dx^2) ;
I = speye(size(L)) ;
 
% Modified elliptic operator
M = E1*L - gamma*(I-E);
  
%% Deal with the boundary
% Compute the modified closest points which will be used to deal with boundary
[cpbarx, cpbary] = cpbar_2d(x, y, cpf{i});

% Build the interpolating matrix for ghost points
E_bdy=interp2_matrix(x1d,y1d,cpbarx(bdy), cpbary(bdy), p,band);

% Modidy the entries of the operator matrix
M(bdy,:)=(I(bdy,:)+E_bdy)/2;

%% Build RHS
rhs = f(cpx, cpy);
    
% bilinear interpolation from the grid points to the closest points of the
% ghost_points 
E_ghost = interp2_matrix(x1d, y1d,cpx(bdy), cpy(bdy), 1, band);

% Impose the boundary condition
rhs(bdy) = E_ghost*g(x,y); 

d{i}.band = band;
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
end


%% Initial guess u0
u0 = {};
for i=1:N
    s= d{i}; 
    u0{i} = cos(10*s.x) + sin(12*s.y);
end

% 2D plot of the initial guess
figure(N+1); clf;

for i=1:N
    s= d{i}; 
    scatter(s.cpx, s.cpy, 40, u0{i}, 'filled');
    xlim([-1.2, 1.2]);
    ylim([-1.2, 1.2]);
    hold on;
    axis equal
    title('Plot of u^{0}')
    xlabel('x'); ylabel('y'); 
    drawnow
end




%% Find the ghost points for schwarz bc
  for i=1:N
      s = d{i};
      j = setdiff(randperm(N), i);
      o = d{j};
 
    % now "s" is the current object we're solving on and "o" is the
    % other object.  Locate ghost points for s that are in the
    % interior of o.
    [C, Io, Is] = intersect(o.band(~ o.bdy), s.band(s.bdy));  % find the ghost points
 
    temp= find(~o.bdy);  Io = temp(Io);
    temp = find(s.bdy);  Is = temp(Is);
    assert (all (C == o.band(Io)));
    assert (all (C == s.band(Is)));
    
     % Plot of the ghost points
     figure(i); 
     scatter(s.x(Is),s.y(Is), 40, 'r')
     axis equal; hold on;   
   
     E_bdy_in=interp2_matrix(x1d, y1d, s.cpbarx(Is), s.cpbary(Is),  1, s.band);
     s.M(Is,:)=(s.I(Is,:)+E_bdy_in)/2;
% Trilinear interpolation from the grid points to the closest points of the
% ghost_points (with respect to the other sub-domains)

     E_ghost_in= interp2_matrix(x1d, y1d,  s.cpx(Is), s.cpy(Is), 1, o.band);
     
     d{i}.Is=Is;
     d{i}.E_ghost_in=E_ghost_in;
  end
    
 %% Schwarz iterations

Nsteps=20; %Number of iterations
 for step=1:Nsteps;
      
     if step == 1
      %d{1}.u=u0{1};  %parallel computing 
       d{2}.u=u0{2};
     end
      
 % Update schwarz b.c. on subdomain 1
 d{1}.rhs(d{1}.Is) = d{1}.E_ghost_in*d{2}.u;
  
 % Solve on subdomain 1
 tic
 d{1}.u= (d{1}.M)\d{1}.rhs;
 toc
 
 % Update schwarz b.c. on subdomain 1 
 d{2}.rhs(d{2}.Is)=d{2}.E_ghost_in*d{1}.u;
 
 % Solve on subdomain 1
 tic
 d{2}.u = (d{2}.M)\d{2}.rhs;
 toc
 
 
 
% Parallel computing: additive schwarz 
 
% d{1}.rhs(d{1}.Is) = d{1}.E_ghost_in*d{2}.u;
% d{2}.rhs(d{2}.Is)=d{2}.E_ghost_in*d{1}.u;
% d{1}.u= (d{1}.M)\d{1}.rhs;
% d{2}.u = (d{2}.M)\d{2}.rhs;

 err = [step];
 
 
  
 %TODO: recompute the error based on the original algorithm
 for i=1:N;
     s=d{i};
  ue_in = uexact(s.x(~(s.bdy)), s.y(~(s.bdy)));
  u_in = s.u(~(s.bdy));
 err(i+1) = max(abs(ue_in-u_in));
 relerr=max(abs(ue_in-u_in))/(max(abs(ue_in)))
 end
 err
 
 
 end
% 
%  
% uexact=(1/4)* (x.^2+y.^2-1);
% f=1;
% g=0;
% overlapping region from x=-0.1 to x=0.1

% TODO: complete the following table
% 
% dx         Max. error     Relative error      Number of iterations
% ------------------------------------------------------------------------------------------
% 0.4       0.0054             0.0365           8
% 0.2        0.0021             0.0128           6
% 0.1        6.4370e-04         0.0039           10
%
