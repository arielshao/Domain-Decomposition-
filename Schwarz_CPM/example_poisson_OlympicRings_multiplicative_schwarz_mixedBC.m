% This code aims to implement the Multiplicative Schwarz Method to 
% solve the Poisson equation -\Delta u=f on five mutually overlapping rings
% with mixed boundary conditions.

% Outer dirichlet boundary condition : u=g
% Inner neumann boundary condition : \frac{\partial u} {\partial n} =h


% This code deals with the outer dirichlet boundary condtions using 
% the extension approach, that is,
% u(x_g)=2u(cp(x_g))-u(cpbar(x_g)) where cpbar=cp(2cp(x)-x)
% while imposing the interior dirichlet boundary conditions using the 
% direct approach, that is,
% u_left(ghost_left)=u_right(ghost_left) and
% u_right(ghost_right)=u_left(ghost_right).

% For the inner neumann boundary condition, we do the following:
% u(x_g)=u(cpbar(x_g)))+2h(cp(x_g))
%       = u(cpbar(x_g))+2 (u_x*n_1(cp(x_g))+u_y*n_2(cp(x_g)))




% Include the cp_matrices folder (edit as appropriate)
addpath('../cp_matrices');
% add functions for finding the closest points
addpath('../surfaces');

%% Number of Iterations
t=cputime;
tic;
Nsteps=50; % Number of iterations
T=100*eps; % Time delay
% Generate the vector for errors
errvals_inf=zeros(1, Nsteps);

%% Plot the composite domain 
R=1;
r=0.85;
LW='LineWidth';
% Plot the centre disc
thetas=linspace(0,2*pi,100)';
r_out=R*ones(size(thetas));
r_in=r*ones(size(thetas));
[x_disc_centre_out,y_disc_centre_out]=pol2cart(thetas,r_out);
[x_disc_centre_in,y_disc_centre_in]=pol2cart(thetas,r_in);
x_disc_centre_out=x_disc_centre_out(:); y_disc_centre_out=y_disc_centre_out(:);
x_disc_centre_in=x_disc_centre_in(:); y_disc_centre_in=y_disc_centre_in(:);
figure(1);clf;
plot(x_disc_centre_out,y_disc_centre_out, 'k-', LW,2);
hold on;
plot(x_disc_centre_in,y_disc_centre_in, 'k-', LW,2);
hold on;

x_disc_left1_out=x_disc_centre_out-1.5;
y_disc_left1_out=y_disc_centre_out-0.5;
x_disc_left1_in=x_disc_centre_in-1.5;
y_disc_left1_in=y_disc_centre_in-0.5;
plot(x_disc_left1_out,y_disc_left1_out, 'k-', LW,2);
plot(x_disc_left1_in,y_disc_left1_in, 'k-', LW,2);
hold on;

x_disc_left2_out=x_disc_centre_out-3;
y_disc_left2_out=y_disc_centre_out;
x_disc_left2_in=x_disc_centre_in-3;
y_disc_left2_in=y_disc_centre_in;
plot(x_disc_left2_out,y_disc_left2_out, 'k-', LW,2);
plot(x_disc_left2_in,y_disc_left2_in, 'k-', LW,2);
hold on;

x_disc_right1_out=x_disc_centre_out+1.5;
y_disc_right1_out=y_disc_centre_out-0.5;
x_disc_right1_in=x_disc_centre_in+1.5;
y_disc_right1_in=y_disc_centre_in-0.5;
plot(x_disc_right1_out,y_disc_right1_out, 'k-', LW,2);
plot(x_disc_right1_in,y_disc_right1_in, 'k-', LW,2);
hold on;

x_disc_right2_out=x_disc_centre_out+3;
y_disc_right2_out=y_disc_centre_out;
x_disc_right2_in=x_disc_centre_in+3;
y_disc_right2_in=y_disc_centre_in;
plot(x_disc_right2_out,y_disc_right2_out, 'k-', LW,2);
plot(x_disc_right2_in,y_disc_right2_in, 'k-', LW,2);
hold on;
xlim([-5, 5])
ylim([-2, 2])


 %% Generate grid points around the composite domain 
dx=0.05; % grid size
x1d=(-4-5*dx:dx:4+5*dx)';
y1d=(-1.5-5*dx:dx:1+5*dx)';
% Generate meshgrid
[xx,yy]=meshgrid(x1d,y1d);
% Make into vectors
X=xx(:); Y=yy(:);
% Plot the grid points
plot(X, Y, '.', 'color', 0.75*[1 1 1]);
hold on;

%% Closest Point representations
% CP points with respect to the rings 1,2,3,4,5 (numbered from left to right)
[cp_x1, cp_y1,dist1,bdy1]=cpAnnulus(X,Y,r,R,[-3,0]);
plot(cp_x1,cp_y1, 'bo')   
[cp_x2, cp_y2,dist2,bdy2]=cpAnnulus(X,Y,r,R,[-1.5,-0.5]);
plot(cp_x2,cp_y2, 'yo')
[cp_x3, cp_y3,dist3,bdy3]=cpAnnulus(X,Y,r,R,[0,0]);
plot(cp_x3,cp_y3, 'ko')
[cp_x4, cp_y4,dist4,bdy4]=cpAnnulus(X,Y,r,R,[1.5,-0.5]);
plot(cp_x4,cp_y4, 'go')
[cp_x5, cp_y5,dist5,bdy5]=cpAnnulus(X,Y,r,R,[3,0]);
plot(cp_x5,cp_y5, 'ro')

%Find the modified closest point representation of each grid point
% "cpbar" [Macdonald, Brandman, Ruuth 2011]: cpbar(x):=cp(2*cp(x)-x);
[cpbar_x1, cpbar_y1, dist6]=cpAnnulus(2*cp_x1-X, 2*cp_y1-Y,...
    r,R,[-3,0]);
plot(cpbar_x1,cpbar_y1,'bo')
[cpbar_x2, cpbar_y2, dist7]=cpAnnulus(2*cp_x2-X, 2*cp_y2-Y,...
    r,R,[-1.5,-0.5]);
plot(cpbar_x2,cpbar_y2,'yo')
[cpbar_x3, cpbar_y3, dist8]=cpAnnulus(2*cp_x3-X, 2*cp_y3-Y,...
   r,R,[0,0]);
 plot(cpbar_x3,cpbar_y3,'ko')
[cpbar_x4, cpbar_y4, dist9]=cpAnnulus(2*cp_x4-X, 2*cp_y4-Y,...
    r,R,[1.5,-0.5]);
plot(cpbar_x4,cpbar_y4,'go')
[cpbar_x5, cpbar_y5, dist10]=cpAnnulus(2*cp_x5-X, 2*cp_y5-Y,...
    r,R,[3, 0]);
plot(cpbar_x5,cpbar_y5,'ro')
%% Initial Condition
% initguess=@(x,y) abs(log(1./((x.^2+y.^2).^.5))).^0.4999; 
% initguess=@(x,y) ones(size(x));
% initguess = @(x,y) zeros(size(x)); % u
% initguess = @(x,y) -1-sin(22*x).*sin(10*y);
 initguess = @(x,y) 1-sin(22*x).*sin(10*y);
%  initguess = @(x,y) 0.0+(min(x.^2+(y+.5).^2,(x-.5).^2+y.^2)<.005);
%  initguess = @(x,y) (min(x.^2+(y+.5).^2,(x-.5).^2+y.^2)>.05);
%  initguess = @(x,y) -3+4*(x.^2+(y+1/2).^2>.001);
%  initguess = @(x,y) 1-.001./(x.^2+(y+1/2).^2).^(1/3)-.001./((x-1/2).^2+y.^2).^(1/3);
%  initguess = @(x,y) 10-7.*(abs(log(-x.^2-(y+1/2).^2))).^0.4999;


% This makes the initial guess u0 into a vector
 u0 = initguess(X,Y);
 %% Assign exact solution 

 uexact= @(x,y)(1/4)*(1- x.^2-y.^2);
 uxexact=@(x,y)-x/2;
 uyexact=@(x,y)-y/2;
 
% uexact=@(x,y) 1+0*x;
u1_exact =uexact(cp_x1,cp_y1);
u2_exact=uexact(cp_x2,cp_y2);
u3_exact =uexact(cp_x3,cp_y3);
u4_exact=uexact(cp_x4,cp_y4);
u5_exact =uexact(cp_x5,cp_y5);


 %% f values
% syms x y;
% u=((x+3).^2+y.^2-1).*((x+1.5).^2+(y+0.5).^2-1).*(x.^2+y.^2-1).*...
%     ((x-1.5).^2+(y+0.5).^2-1).*((x-3).^2+y.^2-1);
% f=-diff(u,x,2)-diff(u,y,2);
% f=matlabFunction(f);
% f_1=f(X,Y);
% f_2=f(X,Y);
% f_3=f(X,Y);
% f_4=f(X,Y);
% f_5=f(X,Y);
% f_1=zeros(size(X));
% f_2=zeros(size(X));
% f_3=zeros(size(X));
% f_4=zeros(size(X));
% f_5=zeros(size(X));
f_1=ones(size(X));
f_2=ones(size(X));
f_3=ones(size(X));
f_4=ones(size(X));
f_5=ones(size(X));
%% Inner boundary condtions

% Neumann boundary conditions
% u_Gamma_1_in=zeros(size(X)); 
% u_Gamma_2_in=zeros(size(X)); 
% u_Gamma_3_in=zeros(size(X)); 
% u_Gamma_4_in=zeros(size(X)); 
% u_Gamma_5_in=zeros(size(X)); 
% u_Gamma_1_in=ones(size(X)); 
% u_Gamma_2_in=ones(size(X)); 
% u_Gamma_3_in=ones(size(X)); 
% u_Gamma_4_in=ones(size(X)); 
% u_Gamma_5_in=ones(size(X)); 

nx_1=X-cp_x1;ny_1=Y-cp_y1;
nx_2=X-cp_x2;ny_2=Y-cp_y2;
nx_3=X-cp_x3;ny_3=Y-cp_y3;
nx_4=X-cp_x4;ny_4=Y-cp_y4;
nx_5=X-cp_x5;ny_5=Y-cp_y5;

u_Gamma_1_in=nx_1.*uxexact(cp_x1,cp_y1)+ny_1.*uyexact(cp_x1,cp_y1);
u_Gamma_2_in=nx_2.*uxexact(cp_x2,cp_y2)+ny_2.*uyexact(cp_x2,cp_y2);
u_Gamma_3_in=nx_3.*uxexact(cp_x3,cp_y3)+ny_3.*uyexact(cp_x3,cp_y3);
u_Gamma_4_in=nx_4.*uxexact(cp_x4,cp_y4)+ny_4.*uyexact(cp_x4,cp_y4);
u_Gamma_5_in=nx_5.*uxexact(cp_x5,cp_y5)+ny_5.*uyexact(cp_x5,cp_y5);


%% Outer boundary conditions
% u_Gamma_1_out=zeros(size(X));
% u_Gamma_2_out=zeros(size(X));
% u_Gamma_3_out=zeros(size(X));
% u_Gamma_4_out=zeros(size(X));
% u_Gamma_5_out=zeros(size(X));

% u_Gamma_1_out=ones(size(X));
% u_Gamma_2_out=ones(size(X));
% u_Gamma_3_out=ones(size(X));
% u_Gamma_4_out=ones(size(X));
% u_Gamma_5_out=ones(size(X));

% Assign the 
u_Gamma_1_out=u1_exact;
u_Gamma_2_out=u2_exact;
u_Gamma_3_out=u3_exact;
u_Gamma_4_out=u4_exact;
u_Gamma_5_out=u5_exact;


%% Build Laplacian and Interpolation matrices
%Construct an interpolation matrix for closest points
dim= 2; % dimension
p=3; % interpolation degree
disp('building laplacian and interp matrices');
[Dxx,Dyy,Dxc,Dyc,Dxb,Dyb,Dxf,Dyf,Dxyc]=diff2d_matrices(x1d,y1d);
L=Dxx+Dyy;
E_1=interp2_matrix(x1d,y1d, cp_x1, cp_y1, p);
E_2=interp2_matrix(x1d,y1d, cp_x2, cp_y2, p);
E_3=interp2_matrix(x1d,y1d, cp_x3, cp_y3, p);
E_4=interp2_matrix(x1d,y1d, cp_x4, cp_y4, p);
E_5=interp2_matrix(x1d,y1d, cp_x5, cp_y5, p);
E1_1=interp2_matrix(x1d,y1d, cp_x1, cp_y1, 1);
E1_2=interp2_matrix(x1d,y1d, cp_x2, cp_y2, 1);
E1_3=interp2_matrix(x1d,y1d, cp_x3, cp_y3, 1);
E1_4=interp2_matrix(x1d,y1d, cp_x4, cp_y4, 1);
E1_5=interp2_matrix(x1d,y1d, cp_x5, cp_y5, 1);
I= speye(size(L));
% Modified Laplacian operator.
% Details can be found in [Macdonald & Brandman & Ruuth 2011]
L_1= E1_1*L-2*dim/dx^2*(I-E_1);
L_2= E1_2*L-2*dim/dx^2*(I-E_2);
L_3= E1_3*L-2*dim/dx^2*(I-E_3);
L_4= E1_4*L-2*dim/dx^2*(I-E_4);
L_5= E1_5*L-2*dim/dx^2*(I-E_5);
M_1=-L_1;
M_2=-L_2;
M_3=-L_3;
M_4=-L_4;
M_5=-L_5;
disp('done');
 
 %% Find the ghost points we need for the Dirichlet boundary conditions
 
%ring1, ring2, ring3, ring4, and ring5 give the set of grid points within
%each ring respectively


tol=dx/10;
ring1=(sqrt((cpbar_x1-X).^2+(cpbar_y1-Y).^2)-sqrt((X-cp_x1).^2+(Y-cp_y1).^2)<tol);
ring2=(sqrt((cpbar_x2-X).^2+(cpbar_y2-Y).^2)-sqrt((X-cp_x2).^2+(Y-cp_y2).^2)<tol);
ring3=(sqrt((cpbar_x3-X).^2+(cpbar_y3-Y).^2)-sqrt((X-cp_x3).^2+(Y-cp_y3).^2)<tol);  
ring4=(sqrt((cpbar_x4-X).^2+(cpbar_y4-Y).^2)-sqrt((X-cp_x4).^2+(Y-cp_y4).^2)<tol);
ring5=(sqrt((cpbar_x5-X).^2+(cpbar_y5-Y).^2)-sqrt((X-cp_x5).^2+(Y-cp_y5).^2)<tol);

% plot(X(ring1),Y(ring1),'ro')
%ringNumber_band (Number=1,2,3,4,5) give the set of grid points within the
%neigbourhood of each ring respectively
% ring1_band=(dist1 <=2*dx);
% ring2_band=(dist2 <=2*dx);
% ring3_band=(dist3 <=2*dx);
% ring4_band=(dist4 <=2*dx);
% ring5_band=(dist5 <=2*dx);

%ringNumber_band_out give the set of grid points on ther outer
%neighbourhood of each ring;
%ringNumber_band_in give the set of grid points on ther inner
%neighbourhood of each ring;
ring1_band_out=((X+3).^2+Y.^2<=R.^2+3*dx & (X+3).^2+Y.^2>=r.^2);
ring2_band_out=((X+1.5).^2+(Y+0.5).^2<=R.^2+3*dx & (X+1.5).^2+(Y+0.5).^2>=r.^2);
ring3_band_out=(X.^2+Y.^2<=R.^2+3*dx & X.^2+Y.^2>=r.^2);
ring4_band_out=((X-1.5).^2+(Y+0.5).^2<=R.^2+3*dx & (X-1.5).^2+(Y+0.5).^2>=r.^2);
ring5_band_out=((X-3).^2+Y.^2<=R.^2+3*dx & (X-3).^2+Y.^2>=r.^2);


ring1_band_in=((X+3).^2+Y.^2<=R.^2 & (X+3).^2+Y.^2>=r.^2-3*dx);
ring2_band_in=((X+1.5).^2+(Y+0.5).^2<=R.^2 & (X+1.5).^2+(Y+0.5).^2>=r.^2-3*dx);
ring3_band_in=(X.^2+Y.^2<=R.^2 & X.^2+Y.^2>=r.^2-3*dx);
ring4_band_in=((X-1.5).^2+(Y+0.5).^2<=R.^2 & (X-1.5).^2+(Y+0.5).^2>=r.^2-3*dx);
ring5_band_in=((X-3).^2+Y.^2<=R.^2 & (X-3).^2+Y.^2>=r.^2-3*dx);


%ghost_Gamma_Number_out gives the set of points we need to impose outer dirichlet
%boundary

ghost_Gamma_1_out=(ring1_band_out & ~ring2 & ~ring1);
ghost_Gamma_2_out=(ring2_band_out & ~ring1 & ~ring2 &~ring3);
ghost_Gamma_3_out=(ring3_band_out & ~ring2 & ~ring3 &~ring4);
ghost_Gamma_4_out=(ring4_band_out & ~ring3 & ~ring4 &~ring5);
ghost_Gamma_5_out=(ring5_band_out & ~ring4 & ~ring5);

%ghost_Gamma_Number_out gives the set of points we need to impose inner
%neumann boundary

ghost_Gamma_1_in=(ring1_band_in & ~ring1& ~ring2);
ghost_Gamma_2_in=(ring2_band_in & ~ring1& ~ring2 &~ring3);
ghost_Gamma_3_in=(ring3_band_in & ~ring2& ~ring3 & ~ring4);
ghost_Gamma_4_in=(ring4_band_in & ~ring3& ~ring4 & ~ring5);
ghost_Gamma_5_in=(ring5_band_in & ~ring4& ~ring5);


ghost_12=(~ring1 & ring2 & (ring1_band_out|ring1_band_in));
ghost_21=(~ring2 & ring1 & (ring2_band_out|ring2_band_in));
ghost_23=(~ring2 & ring3 & (ring2_band_out|ring2_band_in));
ghost_34=(~ring3 & ring4 & (ring3_band_out|ring3_band_in));
ghost_32=(~ring3 & ring2 & (ring3_band_out|ring3_band_in));
ghost_43=(~ring4 & ring3 & (ring4_band_out|ring4_band_in));
ghost_45=(~ring4 & ring5 & (ring4_band_out|ring4_band_in));
ghost_54=(~ring5 & ring4 & (ring5_band_out|ring5_band_in));


%% Building matrices used to impose dirichlet boundary conditions
disp('building matrices to impose dirichlet boundary conditions')

%Interpolation matrices 
E12_bar=interp2_matrix(x1d,y1d,cpbar_x1(ghost_12),cpbar_y1(ghost_12),p);
E21_bar=interp2_matrix(x1d,y1d,cpbar_x2(ghost_21),cpbar_y2(ghost_21),p);
E23_bar=interp2_matrix(x1d,y1d,cpbar_x2(ghost_23),cpbar_y2(ghost_23),p);
E32_bar=interp2_matrix(x1d,y1d,cpbar_x3(ghost_32),cpbar_y3(ghost_32),p);
E34_bar=interp2_matrix(x1d,y1d,cpbar_x3(ghost_34),cpbar_y3(ghost_34),p);
E43_bar=interp2_matrix(x1d,y1d,cpbar_x4(ghost_43),cpbar_y4(ghost_43),p);
E45_bar=interp2_matrix(x1d,y1d,cpbar_x4(ghost_45),cpbar_y4(ghost_45),p);
E54_bar=interp2_matrix(x1d,y1d,cpbar_x5(ghost_54),cpbar_y5(ghost_54),p);
E_Gamma_1_out_bar=interp2_matrix(x1d,y1d,cpbar_x1(ghost_Gamma_1_out),cpbar_y1(ghost_Gamma_1_out),p);
E_Gamma_2_out_bar=interp2_matrix(x1d,y1d,cpbar_x2(ghost_Gamma_2_out),cpbar_y2(ghost_Gamma_2_out),p);
E_Gamma_3_out_bar=interp2_matrix(x1d,y1d,cpbar_x3(ghost_Gamma_3_out),cpbar_y3(ghost_Gamma_3_out),p);
E_Gamma_4_out_bar=interp2_matrix(x1d,y1d,cpbar_x4(ghost_Gamma_4_out),cpbar_y4(ghost_Gamma_4_out),p);
E_Gamma_5_out_bar=interp2_matrix(x1d,y1d,cpbar_x5(ghost_Gamma_5_out),cpbar_y5(ghost_Gamma_5_out),p);
E_Gamma_1_in_bar=interp2_matrix(x1d,y1d,cpbar_x1(ghost_Gamma_1_in),cpbar_y1(ghost_Gamma_1_in),p);
E_Gamma_2_in_bar=interp2_matrix(x1d,y1d,cpbar_x2(ghost_Gamma_2_in),cpbar_y2(ghost_Gamma_2_in),p);
E_Gamma_3_in_bar=interp2_matrix(x1d,y1d,cpbar_x3(ghost_Gamma_3_in),cpbar_y3(ghost_Gamma_3_in),p);
E_Gamma_4_in_bar=interp2_matrix(x1d,y1d,cpbar_x4(ghost_Gamma_4_in),cpbar_y4(ghost_Gamma_4_in),p);
E_Gamma_5_in_bar=interp2_matrix(x1d,y1d,cpbar_x5(ghost_Gamma_5_in),cpbar_y5(ghost_Gamma_5_in),p);
 % Change entries of M
%  M_ghost_12=(I(ghost_12,:)+E12_bar)/2;
%  M_ghost_21=(I(ghost_21,:)+E21_bar)/2;
%  M_ghost_23=(I(ghost_23,:)+E23_bar)/2;
%  M_ghost_32=(I(ghost_32,:)+E32_bar)/2;
%  M_ghost_34=(I(ghost_33,:)+E34_bar)/2;
%  M_ghost_43=(I(ghost_43,:)+E43_bar)/2;
%  M_ghost_45=(I(ghost_45,:)+E45_bar)/2;
%  M_ghost_54=(I(ghost_54,:)+E54_bar)/2;
M_ghost_12=I(ghost_12,:);
M_ghost_21=I(ghost_21,:);
M_ghost_23=I(ghost_23,:);
M_ghost_32=I(ghost_32,:);
M_ghost_34=I(ghost_34,:);
M_ghost_43=I(ghost_43,:);
M_ghost_45=I(ghost_45,:);
M_ghost_54=I(ghost_54,:);

M_ghost_Gamma_1_out=(I(ghost_Gamma_1_out,:)+E_Gamma_1_out_bar)/2;
M_ghost_Gamma_2_out=(I(ghost_Gamma_2_out,:)+E_Gamma_2_out_bar)/2;
M_ghost_Gamma_3_out=(I(ghost_Gamma_3_out,:)+E_Gamma_3_out_bar)/2;
M_ghost_Gamma_4_out=(I(ghost_Gamma_4_out,:)+E_Gamma_4_out_bar)/2;
M_ghost_Gamma_5_out=(I(ghost_Gamma_5_out,:)+E_Gamma_5_out_bar)/2;


M_ghost_Gamma_1_in=(I(ghost_Gamma_1_in,:)-E_Gamma_1_in_bar)/2;
M_ghost_Gamma_2_in=(I(ghost_Gamma_2_in,:)-E_Gamma_2_in_bar)/2;
M_ghost_Gamma_3_in=(I(ghost_Gamma_3_in,:)-E_Gamma_3_in_bar)/2;
M_ghost_Gamma_4_in=(I(ghost_Gamma_4_in,:)-E_Gamma_4_in_bar)/2;
M_ghost_Gamma_5_in=(I(ghost_Gamma_5_in,:)-E_Gamma_5_in_bar)/2;
% 
% 
M_1(ghost_12,:)=M_ghost_12;
M_2(ghost_21,:)=M_ghost_21;
M_2(ghost_23,:)=M_ghost_23;
M_3(ghost_32,:)=M_ghost_32;
M_3(ghost_34,:)=M_ghost_34;
M_4(ghost_43,:)=M_ghost_43;
M_4(ghost_45,:)=M_ghost_45;
M_5(ghost_54,:)=M_ghost_54;

M_1(ghost_Gamma_1_in,:)=M_ghost_Gamma_1_in;
M_2(ghost_Gamma_2_in,:)=M_ghost_Gamma_2_in;
M_3(ghost_Gamma_3_in,:)=M_ghost_Gamma_3_in;
M_4(ghost_Gamma_4_in,:)=M_ghost_Gamma_4_in;
M_5(ghost_Gamma_5_in,:)=M_ghost_Gamma_5_in;

M_1(ghost_Gamma_1_out,:)=M_ghost_Gamma_1_out;
M_2(ghost_Gamma_2_out,:)=M_ghost_Gamma_2_out;
M_3(ghost_Gamma_3_out,:)=M_ghost_Gamma_3_out;
M_4(ghost_Gamma_4_out,:)=M_ghost_Gamma_4_out;
M_5(ghost_Gamma_5_out,:)=M_ghost_Gamma_5_out;


% Bilinear interpolation from the grid points to the closest points of the
% % ghost_points (with respect to each disc)
E_ghost_12=interp2_matrix(x1d, y1d, cp_x1(ghost_12),cp_y1(ghost_12),1);
E_ghost_21=interp2_matrix(x1d, y1d, cp_x2(ghost_21),cp_y2(ghost_21),1);
E_ghost_23=interp2_matrix(x1d, y1d, cp_x2(ghost_23),cp_y2(ghost_23),1);
E_ghost_32=interp2_matrix(x1d, y1d, cp_x3(ghost_32),cp_y3(ghost_32),1);
E_ghost_34=interp2_matrix(x1d, y1d, cp_x3(ghost_34),cp_y3(ghost_34),1);
E_ghost_43=interp2_matrix(x1d, y1d, cp_x4(ghost_43),cp_y4(ghost_43),1);
E_ghost_45=interp2_matrix(x1d, y1d, cp_x4(ghost_45),cp_y4(ghost_45),1);
E_ghost_54=interp2_matrix(x1d, y1d, cp_x5(ghost_54),cp_y5(ghost_54),1);

E_ghost_Gamma_1_out=interp2_matrix(x1d, y1d, cp_x1(ghost_Gamma_1_out),cp_y1(ghost_Gamma_1_out),1);
E_ghost_Gamma_2_out=interp2_matrix(x1d, y1d, cp_x2(ghost_Gamma_2_out),cp_y2(ghost_Gamma_2_out),1);
E_ghost_Gamma_3_out=interp2_matrix(x1d, y1d, cp_x3(ghost_Gamma_3_out),cp_y3(ghost_Gamma_3_out),1);
E_ghost_Gamma_4_out=interp2_matrix(x1d, y1d, cp_x4(ghost_Gamma_4_out),cp_y4(ghost_Gamma_4_out),1);
E_ghost_Gamma_5_out=interp2_matrix(x1d, y1d, cp_x5(ghost_Gamma_5_out),cp_y5(ghost_Gamma_5_out),1);

E_ghost_Gamma_1_in=interp2_matrix(x1d, y1d, cp_x1(ghost_Gamma_1_in),cp_y1(ghost_Gamma_1_in),1);
E_ghost_Gamma_2_in=interp2_matrix(x1d, y1d, cp_x2(ghost_Gamma_2_in),cp_y2(ghost_Gamma_2_in),1);
E_ghost_Gamma_3_in=interp2_matrix(x1d, y1d, cp_x3(ghost_Gamma_3_in),cp_y3(ghost_Gamma_3_in),1);
E_ghost_Gamma_4_in=interp2_matrix(x1d, y1d, cp_x4(ghost_Gamma_4_in),cp_y4(ghost_Gamma_4_in),1);
E_ghost_Gamma_5_in=interp2_matrix(x1d, y1d, cp_x5(ghost_Gamma_5_in),cp_y5(ghost_Gamma_5_in),1);

disp ('done');

%% Building the Plotting matrices
disp ('building plotting matrices');

% Grid points on disc 1, 2,3,4 and 5

xx1=X(ring1); yy1=Y(ring1);
% Interpolate from the embedded space to the grid points on the 1st disc 
Eplot1=interp2_matrix(x1d,y1d,xx1,yy1, p);

xx2=X(ring2); yy2=Y(ring2);
% Interpolate from the embedded space to the grid points on the 2nd disc 
Eplot2=interp2_matrix(x1d,y1d,xx2,yy2, p);

xx3=X(ring3); yy3=Y(ring3);
% Interpolate from the embedded space to the grid points on the 3rd disc 
Eplot3=interp2_matrix(x1d,y1d,xx3,yy3, p);

xx4=X(ring4); yy4=Y(ring4);
% Interpolate from the embedded space to the grid points on the 4th disc 
Eplot4=interp2_matrix(x1d,y1d,xx4,yy4, p);

xx5=X(ring5); yy5=Y(ring5);
% Interpolate from the embedded space to the grid points on the 5th disc 
Eplot5=interp2_matrix(x1d,y1d,xx5,yy5, p);


% Grid points on disc 1, 2, 3, 4, 5 but not intersecting its neigbours
index_11=(ring1 & ~ring2);
xxx1=X(index_11); yyy1=Y(index_11);
Eplot_1=interp2_matrix(x1d,y1d,xxx1,yyy1, p);

index_22=(ring2 & ~ring1 & ~ ring3);
xxx2=X(index_22); yyy2=Y(index_22);
Eplot_2=interp2_matrix(x1d,y1d,xxx2,yyy2, p);

index_33=(ring3 & ~ring2 & ~ ring4);
xxx3=X(index_33); yyy3=Y(index_33);
Eplot_3=interp2_matrix(x1d,y1d,xxx3,yyy3, p);

index_44=(ring4 & ~ring3 & ~ ring5);
xxx4=X(index_44); yyy4=Y(index_44);
Eplot_4=interp2_matrix(x1d,y1d,xxx4,yyy4, p);

index_55=(ring5 & ~ring4);
xxx5=X(index_55); yyy5=Y(index_55);
Eplot_5=interp2_matrix(x1d,y1d,xxx5,yyy5, p);
disp ('done');
%% Plot initial solution u0
figure(2);
u0_1=Eplot1*u0;
u0_3=Eplot3*u0;
u0_5=Eplot5*u0;
u0_2=Eplot_2*u0;
u0_4=Eplot_4*u0;
plot2d_compdomain([u0_1;u0_2;u0_3;u0_4;u0_5],[xx1;xxx2;xx3;xxx4;xx5],[yy1;yyy2;yy3;yyy4;yy5],dx,dx,2);
h=colorbar;
 %% Iterations
 counter=0;
 figure(3);
for step =1:Nsteps
    if step == 1
        u2=u0;
        u4=u0;
    end
     
 %Impose the dirichlet boundary conditions
    f_1(ghost_12)=u2(ghost_12);
    f_3(ghost_32)=u2(ghost_32);
    f_3(ghost_34)=u4(ghost_34);
    f_5(ghost_54)=u4(ghost_54);
    
    
    f_1(ghost_Gamma_1_out)=E_ghost_Gamma_1_out*u_Gamma_1_out;
    f_3(ghost_Gamma_3_out)=E_ghost_Gamma_3_out*u_Gamma_3_out;
    f_5(ghost_Gamma_5_out)=E_ghost_Gamma_5_out*u_Gamma_5_out;
    f_1(ghost_Gamma_1_in)=E_ghost_Gamma_1_in*u_Gamma_1_in;
    f_3(ghost_Gamma_3_in)=E_ghost_Gamma_3_in*u_Gamma_3_in;
    f_5(ghost_Gamma_5_in)=E_ghost_Gamma_5_in*u_Gamma_5_in;
 % Solve the pde in disc 1,3 and 5 using backslash
 u1=M_1\f_1;
 u3=M_3\f_3;
 u5=M_5\f_5;
 %% Plot the solution 
% the solution on disc 1, 3 and 5
uplot1=Eplot1*u1;
uplot3=Eplot3*u3;
uplot5=Eplot5*u5;
% extend the solution to the rest of the domain
uplot_2=Eplot_2*u2;
uplot_4=Eplot_4*u4;
% Plot the solution on the entire domain
plot2d_compdomain2([uplot1;uplot_2;uplot3;uplot_4;uplot5],[xx1;xxx2;xx3;xxx4;xx5],[yy1;yyy2;yy3;yyy4;yy5],dx,dx,3);
colorbar;

hold on;
% Plot the 5 rings
plot(x_disc_centre_out,y_disc_centre_out, 'k-', LW,2);
plot(x_disc_centre_in,y_disc_centre_in, 'k-', LW,2);
plot(x_disc_left1_out,y_disc_left1_out, 'k-', LW,2);
plot(x_disc_left1_in,y_disc_left1_in, 'k-', LW,2);
plot(x_disc_left2_out,y_disc_left2_out, 'k-', LW,2);
plot(x_disc_left2_in,y_disc_left2_in, 'k-', LW,2);
plot(x_disc_right1_out,y_disc_right1_out, 'k-', LW,2);
plot(x_disc_right1_in,y_disc_right1_in, 'k-', LW,2);
plot(x_disc_right2_out,y_disc_right2_out, 'k-', LW,2);
plot(x_disc_right2_in,y_disc_right2_in , 'k-', LW,2);
xlim([-5, 5])
ylim([-2, 2])
pause(T);

% Update the dirichlet boundary condition

f_2(ghost_21)=u1(ghost_21);
f_2(ghost_23)=u3(ghost_23);
f_4(ghost_43)=u3(ghost_43);
f_4(ghost_45)=u5(ghost_45);
f_2(ghost_Gamma_2_out)=E_ghost_Gamma_2_out*u_Gamma_2_out;
f_4(ghost_Gamma_4_out)=E_ghost_Gamma_4_out*u_Gamma_4_out;
f_2(ghost_Gamma_2_in)=E_ghost_Gamma_2_in*u_Gamma_2_in;
f_4(ghost_Gamma_4_in)=E_ghost_Gamma_4_in*u_Gamma_4_in;
u2=M_2\f_2;
u4=M_4\f_4;

% plot the solution on disc 2, 4 
uplot2=Eplot2*u2;
uplot4=Eplot4*u4;
% extend the solution to the rest of the domain
uplot_1=Eplot_1*u1;
uplot_3=Eplot_3*u3;
uplot_5=Eplot_5*u5;
% Plot the solution on the entire domain
plot2d_compdomain2([uplot_1;uplot2;uplot_3;uplot4;uplot_5],[xxx1;xx2;xxx3;xx4;xxx5],[yyy1;yy2;yyy3;yy4;yyy5],dx,dx,3);
colorbar;

hold on;
% Plot the 5 rings
plot(x_disc_centre_out,y_disc_centre_out, 'k-', LW,2);
plot(x_disc_centre_in,y_disc_centre_in, 'k-', LW,2);
plot(x_disc_left1_out,y_disc_left1_out, 'k-', LW,2);
plot(x_disc_left1_in,y_disc_left1_in, 'k-', LW,2);
plot(x_disc_left2_out,y_disc_left2_out, 'k-', LW,2);
plot(x_disc_left2_in,y_disc_left2_in, 'k-', LW,2);
plot(x_disc_right1_out,y_disc_right1_out, 'k-', LW,2);
plot(x_disc_right1_in,y_disc_right1_in, 'k-', LW,2);
plot(x_disc_right2_out,y_disc_right2_out, 'k-', LW,2);
plot(x_disc_right2_in,y_disc_right2_in , 'k-', LW,2);
xlim([-5, 5]);
ylim([-2, 2]);
pause(T);
counter=counter+1
%% The error calculated with respect to infinity norm
 errvals_inf(step)=max(max(abs([Eplot1*u1; Eplot_2*u2;Eplot3*u3;Eplot_4*u4;Eplot5*u5]-...
         [Eplot1*u1_exact; Eplot_2*u2_exact;Eplot3*u3_exact;Eplot_4*u4_exact;Eplot5*u5_exact])));
  end
 figure(4);clf(4);
nn=1:Nsteps;
semilogy(nn,errvals_inf);
xlabel('Number of Iterations'); % x-axis label
ylabel('errvals_{inf}'); % y-axis label
title ('Plot of errvals_{inf} vs number of iterations');
toc
e=cputime-t;

