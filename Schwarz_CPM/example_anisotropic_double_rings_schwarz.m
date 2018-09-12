%% Solving an anisotropic elliptic equation on a two-ring domain using Schwarz iterations 

%  This example solves an anisotropic elliptic equation on a domain
%  consisting of two overlapping rings with dirichlet boundary condition.

%  The PDE is
%  $$ div [A(y)grad (u(y))] = f$$ in \Omega  
%  $$u=g$$       on\partial\Omega 
%  where \Omega is the composite domain consisting of two overlapping rings

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
% g=uexact;
 %% Build the mesh grid 
dx =0.02;

% make vectors of x, y, positions of the grid
x1d = ((-1-6*dx):dx:(2+6*dx))';
y1d = ((-1.5-6*dx):dx:(1.5+6*dx))';
[xx, yy] = meshgrid(x1d, y1d);
%% Banding
dim = 2;
order = 2;
p = 3;
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));

%% Closest point representations of the 2D subdomains 
cpf = {};
R=1;r=0.8;
%closest points of the grid points
cpf{1}=@(x,y) cpAnnulus(x,y,r,R ,[0 0]);
cpf{2}=@(x,y) cpAnnulus(x,y,r,R, [1 0]);

%% Build the elliptic operator and interpolating matrices for each subdomain

N=length(cpf);
for i=1:N

[cpx, cpy, dist, bdy] = cpf{i}(xx, yy);
band = find(abs(dist) <= bw*dx);  
% % store closest points in the band;
cpx= cpx(band); cpy = cpy(band); 
x = xx(band); y = yy(band); 
bdy = bdy(band);
figure(i);clf;
% Plot the grid points
plot(x, y, '.', 'color', 0.75*[1 1 1]); 
hold on;
axis equal
xlabel('x'); ylabel('y'); 
title (['Closest point representation of the subdomain ' num2str(i)]);
% 2D scatter plot of the cp points
scatter(cpx, cpy, 40, 'filled');
xlim([-2, 2.5]);
ylim([-1.5, 1.5]);
hold on;
% % Build the interpolating matrices
E1 = interp2_matrix(x1d, y1d, cpx, cpy,  1, band);
E = interp2_matrix(x1d, y1d,  cpx, cpy,  p, band);
  
% Build the differential operator and Interpolation matrices
% Generate the function-valued matrix A  
 A={@(x,y)x.^2+1, @(x,y) 0; @(x,y) 0,@(x,y) x.^2+1}
% A={@(x,y)x.^2+1, @(x,y) x; @(x,y) x,@(x,y) y.^2+1}

 % We need our matrix A=[a11(x,y) a12(x,y) ; a21(x,y) a22(x,y)] to be positive-definite
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
% Axf(1,1)=1; Axb(1,1)=1;Ayf(1,1)=1; Ayb(1,1)=1;
% Axf(end,end)=1; Axb(end,end)=1; Ayf(end,end)=1; Ayb(end,end)=1; 

[Dxb,Dxf,Dyb,Dyf]=firstderiv_upw1_2d_matrices(x1d,y1d,band);
% Build the differential operator
L=Dxb*(Axf*P11*Dxf)+Dyb*(Ayf*P22*Dyf);
%  L=Dxb*(Axf*P11*Dxf)+Dyb*(Ayf*P22*Dyf)+Dyb*(Axf*P21*Dxf)+Dxb*(Ayf*P12*Dyf);

%Penalty parameter
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
    
% Bilinear interpolation from the grid points to the closest points of the
% ghost_points 

E_ghost = interp2_matrix(x1d, y1d,cpx(bdy), cpy(bdy), 1, band);

%Impose the boundary condition
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
u0 = {};
for i=1:N
    s= d{i}; 
    u0{i} = cos(10*s.x) + sin(12*s.y);
end

% 2D plot of the initial guess
figure(N+1); clf;

for i=1:N
    s= d{i}; 
    scatter(s.cpx(~s.bdy), s.cpy(~s.bdy), 40, u0{i}(~s.bdy), 'filled');
    xlim([-2, 2.5]);
    ylim([-1.5, 1.5]);
    hold on;
    axis equal
    title('Plot of u^{0}')
    xlabel('x'); ylabel('y'); 
    drawnow
end
hold on;


%% Find the ghost points for schwarz bc
  for i=1:N
      s = d{i};
      j = setdiff(randperm(N), i);
      o = d{j};
 
    % now "s" is the current object we're solving on and "o" is the
    % other object.  Locate ghost points for s that are in the
    % interior of o.
    [C, Io, Is] = intersect(o.band(~o.bdy), s.band(s.bdy));
    temp = find(~o.bdy);  Io = temp(Io);
    temp = find(s.bdy);  Is = temp(Is);
    assert (all (C == o.band(Io)));
    assert (all (C == s.band(Is)));
    
    %Plot of the ghost points
     figure(i+3); 
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
 
 figure(6); clf;
for i=1:N
  s = d{i};
  scatter(s.cpx,s.cpy,40, s.u, 'filled')
  hold on;
  xlim([-2 2.5]);
  ylim([-1.5,1.5]);
  axis equal
  title(['2D plot of u^{' num2str(step) '} for CP Multiplicative Schwarz Method'])
  xlabel('x'); ylabel('y');
  drawnow
end
 end
% 
%  
% uexact= @(x,y)x.^2+y.^2;
% f=@(x,y) 8*x.^2+4;
% g=uexact;   A=[x.^2+1, 0; 0, x.^2+1]
% overlapping region from x=-0.1 to x=0.1

% 
% dx         Max. error     Relative error      Number of iterations
% ------------------------------------------------------------------------------------------
% 0.02        2.38e-3            2.38e-3          9
% 0.01        7.36e-4            7.36e-4          16
% 0.005       2.12e-4            2.12e-4          12


%2nd order accuarate