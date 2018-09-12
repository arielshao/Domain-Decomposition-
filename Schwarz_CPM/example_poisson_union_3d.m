%% Solving a Poisson problem on a composite domain in 3D
% Here we solve
%    $$lap u = f$$
% on a domain consisting of a union of solid 3D objects.
% We use the closest point representation of each object,
% combined with a Schwarz iteration over the objects.
%
% The method of manufactured solutions is used to construct
% an exact solution and thus determine g and h.
% Include the cp_matrices folder (edit as appropriate)
addpath('../cp_matrices');
% add functions for finding the closest points
addpath('../surfaces');
%% Build the mesh grid 
dx = 0.05;

% make vectors of x, y, positions of the grid
x1d = ((-1-6*dx):dx:(1+6*dx))';
y1d = ((-1-6*dx):dx:(1+6*dx))';
z1d = ((-1.5-6*dx):dx:(1.5+6*dx))';
[xx, yy, zz] = meshgrid(x1d, y1d, z1d);

% %% Manufactured solution
%   uexact = @(x,y,z) (exp(sin(4*x)) + x.*sin(5*y))/2 + z;
%   uxexact = @(x,y,z) 2 * exp(sin (4 * x)) .* cos(4 * x) + sin(5 * y) / 2;
%   uyexact = @(x,y,z) 5 * x .* cos(5 * y) / 2;
%   uzexact = @(x,y,z) ones(size(x));
%   f = @(x,y,z) -25 * x .* sin (5 * y) / 2 + 8 * (-sin (4 * x) + cos (4 * x) .^ 2) .* exp (sin (4 * x));

uexact=@(x,y,z) (1/6)* (x.^2+y.^2+z.^2-1);
f=@(x,y,z)0*x+1;
%% Choose the input 3D objects
which = 4;
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
      
    case 6
      cpf{1}=@(x,y,z) cpleftBall(x,y,z, 1, [0 0 0], 0.1);
      cpf{2}=@(x,y,z) cprightBall(x,y,z,1, [0 0 0], -0.1);
end

%% Banding
dim = 3;
order = 2;
p = 3;
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));

%% Loop over each object and build banded grids and operators
N = length(cpf);
for i=1:N
  [cpx, cpy, cpz, dist, bdy] = cpf{i}(xx, yy, zz);
  
  band = find(abs(dist) <= bw*dx);

  % store closest points in the band;
  cpx = cpx(band); cpy = cpy(band); cpz = cpz(band);
  x = xx(band); y = yy(band); z = zz(band);
  [cpbarx, cpbary, cpbarz] = cpbar_3d(x, y, z, cpf{i});

  bdy = bdy(band);


  E1 = interp3_matrix(x1d, y1d, z1d, cpx, cpy, cpz, 1, band);
  E = interp3_matrix(x1d, y1d, z1d, cpbarx, cpbary, cpbarz, p, band);

  E(bdy,:) = -E(bdy,:);
  % note use of non-cpbar
  bc_dir = zeros(size(x));
  bc_dir(bdy) = uexact(cpx(bdy), cpy(bdy), cpz(bdy));
  bc_neu = zeros(size(x));

  %% Create Laplacian matrix
  L = laplacian_3d_matrix(x1d,y1d,z1d,order, band);

  %% Penalty parameter
  gamma = 2*dim/(dx^2) ;
  I = speye(size(L)) ;

  %% Build RHS
  % Note: exact neumann extension done here
  rhs = f(cpx, cpy, cpz);
  %rhs = E1*f(x, y);

  %% This needs to be done later as we adjust the BC
  %rhs = rhs - gamma*(2*bc_dir + 2*bc_neu);

  M = E1*L - gamma*(I-E);

  %% stuff everything into a list of grid "objects"
  g{i}.band = band;
  g{i}.bdy = bdy;
  g{i}.cpx = cpx;
  g{i}.cpy = cpy;
  g{i}.cpz = cpz;
  g{i}.x = x;
  g{i}.y = y;
  g{i}.z = z;
  g{i}.M = M;
  g{i}.bc_dir = bc_dir;
  g{i}.rhs = rhs;
end




%% ICs
uold = {};
for i=1:N
  s = g{i}; % what does this mean?
  uold{i} = cos(10*s.x) + sin(12*s.y) + cos(13*s.z);
end

%% Rough plot
figure(1); clf;
for i=1:N
  s = g{i};
  scatter3(s.cpx, s.cpy, s.cpz, 40, uold{i}, 'filled')
  hold on;
  axis equal
  title('3D plot of u^{0} for CP Multiplicative Schwarz Method')
  xlabel('x'); ylabel('y'); zlabel('z');
  drawnow
end


%% outer iteration (swartz?)
for J = 1:20
u = {};
for i=1:N
  s = g{i};

  %% build a randomly ordered list of the other shapes
  jset = setdiff(randperm(N), i);
  for j=jset
    o = g{j};
    %% now "s" is the current object we're solving on and "o" is the
    % other object.  Locate ghost points for s that are in the
    % interior of o.
    [C, Io, Is] = intersect(o.band(~o.bdy), s.band(s.bdy));
    temp = find(~o.bdy);  Io = temp(Io);
    temp = find(s.bdy);  Is = temp(Is);
    assert (all (C == o.band(Io)));
    assert (all (C == s.band(Is)));
    if (~ isempty(C))
      if (1==1)
        % TODO: bit wasteful to recompute these each time...
	%tic
        Ej = interp3_matrix(x1d, y1d, z1d, s.cpx(Is), s.cpy(Is), s.cpz(Is), 1, o.band);
	%toc
        s.bc_dir(Is) = Ej*uold{j};
      else
        % this was only first-order accurate... Is this the "direct method"?
        s.bc_dir(Is) = uold{j}(Io);
      end
    end
    
    % Plot of the ghost points
%      figure(10+i); clf; 
%      scatter3(s.cpx(Is),s.cpy(Is),s.cpz(Is), 40, 'filled')
%      axis equal; hold on; 
     
    %figure(10+i); clf; I = s.bdy
    %scatter3(xx(I),yy(I),zz(I), 40, zz(I)+20, 'filled')
    %axis equal; hold on; pause
  end

  tic
  rhs = s.rhs - gamma*(2*s.bc_dir);
  u{i} = (s.M)\rhs;
  toc

  %figure(i); clf;
  %I = s.bdy;
  %scatter3(s.cpx(I),s.cpy(I),s.cpz(I), 40, u{i}(I) + 4*(i-1), 'filled')
  %axis equal
  %hold on;
end


err = [J];
for i=1:N
  s = g{i};
  ue_in = uexact(s.x(~(s.bdy)), s.y(~(s.bdy)), s.z(~(s.bdy)));
  u_in = u{i}(~(s.bdy));
  err(i+1) = max(abs(ue_in - u_in));
end
err
uold = u;



% %% slicey dicey plotsies
% [xs1,ys1] = meshgrid(x1d, y1d); zs1 = 0.5*ones(size(xs1));
% [xs2,ys2] = meshgrid(x1d, y1d); zs2 = -0.5*ones(size(xs2));
% [xs3,zs3] = meshgrid(x1d, z1d); ys3 = 0*ones(size(xs3));
% U = 0*xx;
% Ucount = 0*xx;
% for i=1:N
%   figure(10+i); clf;
%   s = g{i};
%   uu = nan*xx;
%   uu(s.band) = u{i};
%   U(s.band(~s.bdy)) = U(s.band(~s.bdy)) + u{i}(~s.bdy);
%   Ucount(s.band(~s.bdy)) = Ucount(s.band(~s.bdy)) + 1;
%   slice(xx,yy,zz,uu,xs1,ys1,zs1)
%   hold on;
%   slice(xx,yy,zz,uu,xs2,ys2,zs2)
%   slice(xx,yy,zz,uu,xs3,ys3,zs3)
%   shading flat
%   xlabel('x'); ylabel('y'); zlabel('z');
%   title(sprintf('subdomain u_{%d}^{(%d)}', i, J))
%   axis equal
%   drawnow
% end
% U = U ./ Ucount;  % average or nan (yuck!)
% figure(10); clf;
% slice(xx,yy,zz,U,xs1,ys1,zs1)
% hold on;
% slice(xx,yy,zz,U,xs2,ys2,zs2)
% slice(xx,yy,zz,U,xs3,ys3,zs3)
% shading flat
% xlabel('x'); ylabel('y'); zlabel('z');
% title(sprintf('average reconstruction u^{(%d)}', J))
% axis equal
% view([-.5 0])
% drawnow



figure(2); clf;
for i=1:N
  s = g{i};
  scatter3(s.cpx,s.cpy, s.cpz, 40, uold{i}, 'filled')
  hold on;
  axis equal
  title(['3D plot of u^{' num2str(J) '} for CP Multiplicative Schwarz Method'])
  xlabel('x'); ylabel('y');zlabel('z');
  drawnow
end
end