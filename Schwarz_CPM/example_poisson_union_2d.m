%% Solving a Poisson problem on a composite domain in 2D
% Here we solve
%    $$lap u = f$$
% on a domain consisting of a union of 2D objects.
% We use the closest point representation of each subdomain,
% combined with a Schwarz iteration over the subdomains.
%
% The method of manufactured solutions is used to construct
% an exact solution and thus determine g and h.
% Include the cp_matrices folder (edit as appropriate)
addpath('../cp_matrices');
% add functions for finding the closest points
addpath('../surfaces');
%% Build the mesh grid 
dx=0.01; % grid size
x1d=(-4-5*dx:dx:4+5*dx)';
y1d=(-1.5-5*dx:dx:1+5*dx)';
% Generate meshgrid
[xx,yy]=meshgrid(x1d,y1d);
% % Make into vectors
% X=xx(:); Y=yy(:);
% % Plot the grid points
% plot(X, Y, '.', 'color', 0.75*[1 1 1]);
% hold on;

 %% Manufactured solution
  uexact = @(x,y) (exp(sin(4*x)) + x.*sin(5*y))/2;
  uxexact = @(x,y) 2 * exp(sin (4 * x)) .* cos(4 * x) + sin(5 * y) / 2;
  uyexact = @(x,y) 5 * x .* cos(5 * y) / 2;
  f = @(x,y) -25 * x .* sin (5 * y) / 2 + 8 * (-sin (4 * x) + cos (4 * x) .^ 2) .* exp (sin (4 * x));
%% Choose the input 2D objects
R=1;
r=0.8; 
which =3;
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
     
  case 5
     cpf{1}=@(x,y)cpleftDisc2(x,y,1,[0,0], 0.2);
     cpf{2}=@(x,y)cprightDisc2(x,y,1,[0,0], -0.2);
end

%% Banding
dim = 2;
order = 2;
p = 3;
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));

%% Loop over each object and build banded grids and operators
N = length(cpf);
for i=1:N
  [cpx, cpy, dist, bdy] = cpf{i}(xx, yy);
  
  band = find(abs(dist) <= bw*dx);

  % store closest points in the band;
  cpx = cpx(band); cpy = cpy(band);
  x = xx(band); y = yy(band);
  [cpbarx, cpbary] = cpbar_2d(x, y,cpf{i});
  bdy = bdy(band);
  E1 = interp2_matrix(x1d, y1d,cpx, cpy,1, band);
  E = interp2_matrix(x1d, y1d,cpbarx, cpbary, p, band);

  E(bdy,:)= -E(bdy,:);
  % note use of non-cpbar
  bc_dir = zeros(size(x));
  bc_dir(bdy) = uexact(cpx(bdy), cpy(bdy));
  bc_neu = zeros(size(x));

  %% Create Laplacian matrix
  L = laplacian_2d_matrix(x1d,y1d,order, band);

  %% Penalty parameter
  gamma = 2*dim/(dx^2) ;
  I = speye(size(L)) ;

  %% Build RHS
  % Note: exact neumann extension done here
  rhs = f(cpx, cpy);
  %rhs = E1*f(x, y);

  %% This needs to be done later as we adjust the BC
  %rhs = rhs - gamma*(2*bc_dir + 2*bc_neu);

  M = E1*L - gamma*(I-E);

  %% stuff everything into a list of grid "objects"
  g{i}.band = band;
  g{i}.bdy = bdy;
  g{i}.cpx = cpx;
  g{i}.cpy = cpy;
  g{i}.x = x;
  g{i}.y = y;
  g{i}.M = M;
  g{i}.bc_dir = bc_dir;
  g{i}.rhs = rhs;
end

%% ICs
uold = {};
for i=1:N
  s = g{i}; % what does this mean?
  uold{i} = cos(10*s.x) + sin(12*s.y);
end

%% Rough plot
figure(1); clf;
for i=1:N
  s = g{i};
  scatter(s.cpx,s.cpy,40, uold{i}, 'filled')
  hold on;
  axis equal
  title('2D plot of u^{0} for CP Multiplicative Schwarz Method')
  xlabel('x'); ylabel('y');
  drawnow
end


%% outer iteration (swartz?)

Nsteps=20;
for J = 1:Nsteps;

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
      if (1==0)
        % TODO: bit wasteful to recompute these each time...
	%tic
        Ej = interp2_matrix(x1d, y1d, s.cpx(Is), s.cpy(Is), 1, o.band);
	%toc
        s.bc_dir(Is) = Ej*uold{j};
      else
        % this was only first-order accurate... Is this the "direct method"?
        s.bc_dir(Is) = uold{j}(Io);
      end
    end
  end

  tic
  rhs = s.rhs - gamma*(2*s.bc_dir);
  u{i} = (s.M)\rhs;
  toc;
   
end

err = [J];
for i=1:N
  s = g{i};
  ue_in = uexact(s.x(~(s.bdy)),s.y(~(s.bdy)));
  u_in = u{i}(~(s.bdy));
  err(i+1) = max(abs(ue_in - u_in));
end

err
uold = u;

figure(2); clf;
for i=1:N
  s = g{i};
  scatter(s.cpx,s.cpy,40, uold{i}, 'filled')
  hold on;
  axis equal
  title(['2D plot of u^{' num2str(J) '} for CP Multiplicative Schwarz Method'])
  xlabel('x'); ylabel('y');
  drawnow
end
end
