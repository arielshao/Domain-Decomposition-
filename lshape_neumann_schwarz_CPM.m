% This code aims to implement the Schwarz Alternating Method to 
% solve a shifted poisson equation -\Delta u+u=f on a L-shape domain 
% without using banding
%The matrices built in this code are adapted from Yujia's code

% Include the cp_matrices folder (edit as appropriate)
addpath('../cp_matrices');
% add functions for finding the closest points
addpath('../surfaces');
%% Number of Iteration
tic;
Nsteps = 80;   % Number of iterations
T = 100*eps;   % Time delay
%% Generate the vector for errors
errvals_inf= zeros(1,Nsteps); 
errvals_H1=zeros(1,Nsteps);
%% Plot the L-shape domain
% Plot the lower-rectangular domain 
xlow1=-1;xlow2=1;
ylow1=-1;ylow2=0;
xlow = [xlow1, xlow2, xlow2, xlow1, xlow1];
ylow = [ylow1, ylow1, ylow2, ylow2, ylow1];
figure(1);clf;
plot(xlow, ylow, 'b-', 'LineWidth', 2);
xlim([-2 2]);
ylim([-2 2]);
hold on;
% Plot the upper-rectangular domain 
xup1=0;xup2=1;
yup1=-1;yup2=1;
xup = [xup1, xup2, xup2, xup1, xup1];
yup = [yup1, yup1, yup2, yup2, yup1];
plot(xup, yup, 'b-', 'LineWidth', 2);
xlim([-2 2]);
ylim([-2 2]);
hold on;

%% Generate grid points around the L-shape domain
dx=0.05; %grid size
x1d = (-1.5:dx:1.5)';
y1d = (-1.5:dx:1.5)';
nx=length(x1d);
ny=length(y1d);
% Generate meshgrid
[x,y]=meshgrid(x1d,y1d);
% Make into vectors
X=x(:);Y=y(:);


%% Initial Condition
initguess=@(x,y) abs(log(1./((x.^2+y.^2).^.5))).^0.4999 % u0=infinity
%  initguess=@(x,y) ones(size(x));
 %initguess = @(x,y) zeros(size(x)); % u0~=1 @(0,0)
% initguess = @(x,y) -1-sin(22*x).*sin(10*y);
%   initguess = @(x,y) 1-sin(22*x).*sin(10*y);
% initguess = @(x,y) 0.0+(min(x.^2+(y+.5).^2,(x-.5).^2+y.^2)<.005);
% initguess = @(x,y) (min(x.^2+(y+.5).^2,(x-.5).^2+y.^2)>.05);
% initguess = @(x,y) -3+4*(x.^2+(y+1/2).^2>.001);
% initguess = @(x,y) 1-.001./(x.^2+(y+1/2).^2).^(1/3)-.001./((x-1/2).^2+y.^2).^(1/3);
% initguess = @(x,y) 10-7.*(abs(log(-x.^2-(y+1/2).^2))).^0.4999;


% This makes the initial guess u0 into a vector
u0 = initguess(X,Y);
%% Closest Point representaitons 
% CP points w.r.t the lower rectangular region
[cpxlow,cpylow,dist1]=cpRectangleDomain(X,Y,2,1,[0,-.5]); 

% CP points w.r.t the upper rectangular region
[cpxup,cpyup,dist2]=cpRectangleDomain(X,Y,1,2,[0.5,0]);

% Find the modified closest point representation of each grid point
% "cpbar" [Macdonald, Brandman, Ruuth 2011]:cpbar(x):=cp(2*cp(x)-x)
[cpbar_xlow, cpbar_ylow,dist3]=cpRectangleDomain(2*cpxlow-X, ...
2*cpylow-Y, 2,1,[0,-.5]); % Lower rectangular region

[cpbar_xup, cpbar_yup,dist4]=cpRectangleDomain(2*cpxup-X, ...
2*cpyup-Y,1,2,[0.5,0]); % Upper rectangular region

%% Build Laplacian and Interpolation matrices

% Construct an interpolation matrix for closest points

dim = 2;  % dimension
p = 3;    % interpolation degree
order=2;  % laplacian order
disp('building laplacian and interp matrices');
[Dxx,Dyy,Dxc,Dyc,Dxb,Dyb,Dxf,Dyf,Dxyc]= diff2d_matrices(x1d,y1d);
Lup=Dxx+Dyy;
Llow=Dxx+Dyy;
Eup = interp2_matrix(x1d,y1d,cpxup, cpyup, p);
Elow=interp2_matrix(x1d,y1d,cpxlow,cpylow,p);
E1up = interp2_matrix(x1d,y1d,cpxup, cpyup, 1);
E1low=interp2_matrix(x1d,y1d,cpxlow,cpylow,1);
Iup = speye(size(Lup));
Ilow=speye(size(Llow));
% Modified Laplacian operator. 
% Details can be found in [Macdonald& Brandman & Ruuth 2011]
Lup= E1up*Iup*Lup-2*dim/dx^2*(Iup - Eup); 
Llow=E1low*Ilow*Llow-2*dim/dx^2*(Ilow-Elow);
Mup=-Lup+Iup;
Mlow=-Llow+Ilow;
disp('done')

%% Find the ghost points we need for the Dirichlet boudnary part
ghost_index=find(sqrt((cpbar_xlow-X).^2+(cpbar_ylow-Y).^2)>...
    sqrt((X-cpxlow).^2+(Y-cpylow).^2) & X>0 &X<1 & Y<1& Y>0)

%ghost1 gives the set of all ghost points
ghost1low=(sqrt((cpbar_xlow-X).^2+(cpbar_ylow-Y).^2)>...
    sqrt((X-cpxlow).^2+(Y-cpylow).^2));
%The following conditions restrict our ghost points to those we need to apply with Dirichlet BCs
ghost2low= (X>=0 & X<=1);
ghost3low= (Y>=0 & Y<=1);
ghost_low=(ghost1low & ghost2low & ghost3low);

%ghost1 gives the set of all ghost points
ghost1up=(sqrt((cpbar_xup-X).^2+(cpbar_yup-Y).^2)>...
    sqrt((X-cpxup).^2+(Y-cpyup).^2));
%The following conditions restrict our ghost points to those we need to apply with Dirichlet BCs
ghost2up= (X<=0 & X>=-1);
ghost3up= (Y<=0 & Y>=-1);
ghost_up=(ghost1up & ghost2up & ghost3up);



%% Building matrices to deal with dirichlet boundary conditions 

disp('buidling matrices to deal with boundary conditions ... ')
Elow_bar = interp2_matrix(x1d,y1d,cpbar_xlow(ghost_low), cpbar_ylow(ghost_low), p);
M_ghost_low = (Ilow(ghost_low,:) + Elow_bar)/2;
Mlow(ghost_low,:) = M_ghost_low;
Eup_bar = interp2_matrix(x1d,y1d,cpbar_xup(ghost_up), cpbar_yup(ghost_up), p);
M_ghost_up = (Iup(ghost_up,:) + Eup_bar)/2;
Mup(ghost_up,:) = M_ghost_up;
disp('done')


%% Building the Plotting matrices
disp('building plotting matrices');
% Grid points on the lower rectangle
x1dlow = (xlow1:dx:xlow2)';
y1dlow= (ylow1:dx:ylow2)';
% Generate the mesh
[xxlow,yylow]=meshgrid(x1dlow,y1dlow);
% Interpolate from the embedded space to the grid points on the lower recntangle only
Eplotlow= interp2_matrix(x1d,y1d,xxlow(:),yylow(:),p);

% Grid points on the upper rectangle
x1dup = (xup1:dx:xup2)';
y1dup= (yup1:dx:yup2)';
% Generate the mesh
[xxup,yyup]=meshgrid(x1dup,y1dup);
% Interpolate from the embedded space to the grid points on the recntangle only
Eplotup= interp2_matrix(x1d,y1d,xxup(:),yyup(:),p);
disp('done');



% Grid points on the upper rectangle but not in the lower one
xx1dup = (xup1:dx:xup2)';
yy1dup = (0:dx:yup2)'; % y>=0
% Generate mesh on the upper square
[xxxup yyyup] = meshgrid(xx1dup, yy1dup);
% Interpolates points to the square
Eplot_up= interp2_matrix(x1d,y1d,xxxup(:),yyyup(:),p);

% Get grid points on the lower rectangle but not in the upper rectangle
xx1dlow = (xlow1:dx:0)';
yy1dlow = (ylow1 :dx:ylow2)'; 
% Generate mesh on the upper rectangle
[xxxlow yyylow] = meshgrid(xx1dlow, yy1dlow);
% Interpolating matrix
Eplot_low= interp2_matrix(x1d,y1d,xxxlow(:),yyylow(:),p);

%% Assign true solution
u_up_true = ones(length(Lup),1); 
u_low_true = ones(length(Llow),1); 
%% Iterations

f_up=ones(length(Lup),1); %RHS
f_low=ones(length(Llow),1)

counter = 0;
for step = 1:Nsteps;
    if step == 1
        uup = u0;
    else
        uup = uup;
    end
% Update the Dirichlet condtion 
f_low(ghost_low)=uup(ghost_low);
 
% Solve the pde in the lower rectangular domain
%(-Lnew+I)u=f for u using backslash.
ulow=Mlow\f_low;
%Plot for the solution resctricted to the lower rectangle
uplotlow=Eplotlow*ulow;
uplotlow=reshape(uplotlow,size(xxlow));
surf(xxlow,yylow, uplotlow);
hold on;

uplot_up=Eplot_up*uup;
uplot_up=reshape(uplot_up,size(xxxup));
% Plot the solution to the rest of the L-shape domain
surf(xxxup,yyyup, uplot_up);
hold on;

pause(T);
counter = counter+1;

clf;

plot(xlow, ylow, 'b-', 'LineWidth', 2);
xlim([-2 2]);
ylim([-2 2]);
hold on;
% Plot the upper-rectangular domain 
plot(xup, yup, 'b-', 'LineWidth', 2);
xlim([-2 2]);
ylim([-2 2]);
hold on;

% Update the Dirichlet condtion vector R_up
f_up(ghost_up)=ulow(ghost_up);
 
% Solve the pde in the up rectangular domain
%(-Lnew+I)u=f for u using backslash.
uup=Mup\f_up;
%Plot for the solution resctricted to the lower rectangle
uplotup=Eplotup*uup;
uplotup=reshape(uplotup,size(xxup));
surf(xxup,yyup, uplotup);
h=colorbar
hold on;


uplot_low=Eplot_low*ulow;
uplot_low=reshape(uplot_low,size(xxxlow));
% Plot the solution to the rest of the L-shape domain
surf(xxxlow,yyylow, uplot_low);
hold on;

pause(T);
counter = counter+1;
% The error calculated with respect to infinity norm
errvals_inf(step)=max(max(abs([uup(:); ulow(:)]-[u_up_true(:); u_low_true(:)])));
end

figure(2)
nn = 1:Nsteps;
semilogy(nn,errvals_inf);
xlabel('Number of Iterations'); % x-axis label
ylabel('errvals_{inf}') ;    % y-axis label
toc





