%% Solving a Poisson's problem on a torus 
%  This example solves the poisson's equation on a 
%  torus

%  The PDE is
%  $$-lap u+c*u = f$$ in \Omega  
%  where \Omega is a sphere and c is a non-negative shift


% The method of manufactured solutions is used to construct
% an exact solution and thus determine f and g.

% adjust as appropriate
addpath('../cp_matrices');
addpath('../surfaces');


R=1;r=0.4;
%% Manufactured solution
c=1;

% uexact=@(th, phi)sin(3*th)+cos(2*phi);
% f=@(th, phi,R,r)-(9*sin(3*th)./((R+r.*cos(phi)).^2)-2*sin(phi).*sin(2*(phi))*r./(R+r.*cos(phi)))...
%     +4*cos(2*phi)./(r.^2))+ c*uexact(th,phi);
 uexact = @(th, phi) sin(3*th) + cos(2*phi);
 f = @(th, phi)  9*sin(3*th)./(R+r*cos(phi)).^2-sin(phi)./(r*(R+r*cos(phi)))*2.*sin(2*phi)+ 4*cos(2*phi)/r^2 +c*uexact(th, phi);
 g=uexact;
%% Build the mesh grid 

dx = 0.1;  %grid size
pad = 5;
x1d=((-3-pad*dx):dx:(3+pad*dx))';
y1d=x1d;
z1d=x1d;

[xx,yy,zz]=meshgrid(x1d,y1d,z1d);


%% Banding
dim = 3;
order = 2;
p = 3;
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));

%% Closest point representations of the 3D subdomains 
cpf = {};
%closest points of the grid points
cpf{1}=@(x,y,z) cpTorus_left(x,y,z,R,r, [0 0 0],  0.1);
cpf{2}=@(x,y,z) cpTorus_right(x,y,z,R, r,[0 0 0], -0.1);

%% Build the elliptic operator and interpolating matrices for each subdomain

N=length(cpf);
for i=1:N
[cpx, cpy, cpz, dist, bdy] = cpf{i}(xx, yy, zz);
band = find(abs(dist) <= bw*dx);
% store closest points in the band;
cpx= cpx(band); cpy = cpy(band); cpz = cpz(band);
x = xx(band); y = yy(band); z = zz(band);
bdy = bdy(band);
[th, phi,r] = cart2paramTorus(x,y,z,R);
[cpth, cpphi,cpr] = cart2paramTorus(cpx,cpy,cpz,R);
  
% 3D scatter plot of the cp points
figure(i);clf;
scatter3(cpx, cpy, cpz, 40, 'filled');
hold on;
axis equal
xlabel('x'); ylabel('y'); zlabel('z');
title (['Closest point representation of subdomain ' num2str(i)]);

% Build the interpolating matrices
E1 = interp3_matrix(x1d, y1d, z1d, cpx, cpy, cpz, 1, band);
E = interp3_matrix(x1d, y1d, z1d, cpx, cpy, cpz, p, band);
  
% Create Laplacian matrix
L = laplacian_3d_matrix(x1d,y1d,z1d,order, band);
 
% Penalty parameter
gamma = 2*dim/(dx^2) ;
I = speye(size(L)) ;
 
% Modified elliptic operator
M = -(E1*L - gamma*(I-E))+c*I;
  
%% Deal with the boundary
% Compute the modified closest points which will be used to deal with boundary
[cpbarx, cpbary, cpbarz] = cpbar_3d(x, y, z, cpf{i});

% Build the interpolating matrix for ghost points
E_bdy=interp3_matrix(x1d,y1d,z1d,cpbarx(bdy), cpbary(bdy),cpbarz(bdy), p,band);

% Modidy the entries of the operator matrix
M(bdy,:)=(I(bdy,:)+E_bdy)/2;

%% Build RHS
% rhs = f(cpth, cpphi, 1);
rhs=f(cpth,cpphi);
    
% Trilinear interpolation from the grid points to the closest points of the
% ghost_points 
E_ghost = interp3_matrix(x1d, y1d, z1d, cpx(bdy), cpy(bdy), cpz(bdy), 1, band);

% Impose the boundary condition
rhs(bdy) = E_ghost*g(th, phi); 

d{i}.band = band;
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
d{i}.th=th;
d{i}.phi=phi;
end


%% Initial guess u0
u0 = {};
for i=1:N
    s= d{i}; 
    u0{i} = cos(10*s.x) + sin(12*s.y) + cos(13*s.z); 
end

% 3D plot of the initial guess
figure(N+1); clf;

for i=1:N
    s= d{i}; 
    scatter3(s.cpx, s.cpy, s.cpz, 40, u0{i}, 'filled');
    hold on;
    axis equal
    title('3D plot of u^{0}')
    xlabel('x'); ylabel('y'); zlabel('z');
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
     figure(N+1+i); clf; 
     scatter3(s.x(Is),s.y(Is),s.z(Is), 40, 'filled')
     axis equal; hold on;   
   
     E_bdy_in=interp3_matrix(x1d, y1d, z1d, s.cpbarx(Is), s.cpbary(Is), s.cpbarz(Is), 1, s.band);
     s.M(Is,:)=(s.I(Is,:)+E_bdy_in)/2;
% Trilinear interpolation from the grid points to the closest points of the
% ghost_points (with respect to the other sub-domains)

     E_ghost_in= interp3_matrix(x1d, y1d, z1d, s.cpx(Is), s.cpy(Is), s.cpz(Is), 1, o.band);
     
     d{i}.Is=Is;
     d{i}.E_ghost_in=E_ghost_in;
  end
    
 %% Schwarz iterations

Nsteps=40; %Number of iterations
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
     ue_in = uexact(s.th(~(s.bdy)), s.phi(~(s.bdy)));
     u_in =  s.u(~(s.bdy));
     err(i+1) = max(abs(ue_in-u_in));
     relerr=max(abs(ue_in-u_in))/(max(abs(ue_in)))
 end
 err
 

 end
%% Plot of the solution after 40 iterations

for i=1:N
    figure(5+i);clf;
    s= d{i}; 
    scatter3(s.cpx(~s.bdy), s.cpy(~s.bdy), s.cpz(~s.bdy), 40, d{i}.u(~s.bdy), 'filled');
    hold on;
    axis equal
    title(['3D scatter plot of u_{' num2str(i) '}^{20}'])
    xlabel('x'); ylabel('y'); zlabel('z');
    drawnow
end

figure(8);clf;

for i=1:N
    s= d{i}; 
    scatter3(s.cpx(~s.bdy), s.cpy(~s.bdy), s.cpz(~s.bdy), 40, d{i}.u(~s.bdy), 'filled');
    hold on;
    axis equal
    title('3D scatter plot of u^{40}')
    xlabel('x'); ylabel('y'); zlabel('z');
    drawnow
    view(-45,45)
end






 

% overlapping region from x=-0.1 to x=0.1

% dx         Error       Relative error    Schwarz: Error     Relative  error     # iterations 
% --------------------------------------------------------------------------------------------        
% 0.2       1.78e-1        8.94e-2             2.05e-1            1.03e-1                21
% 0.1       7.99e-2        4.00e-2             1.27e-1            6.33e-2                7
% 0.05      1.20e-2        6.00e-3             3.20e-2            1.60e-2               29
% 
 
% Special case dx=0.025 : it is probably due to the fact that dx is too
% small, and thus the overlapping region is too large (banding results in
% inaccuaracy) 

% It is 2nd order accuarate 

% TODO : the ghost points are asymmetric (thus problematic),so check the
% ghost points, cpTorus_left; cpTorus_right!!! 
 
