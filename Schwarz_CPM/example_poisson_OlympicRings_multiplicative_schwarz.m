% This code aims to implement the Multiplicative Schwarz Method to 
% solve the Poisson equation -\Delta u=f on five mutually overlapping rings
% with homogeneous or non-homoegenous dirichlet boundary conditions.

% This code deals with the outer dirichlet boundary condtions using 
% the extension approach, that is,
% u(x_g)=2u(cp(x_g))-u(cpbar(x_g)) where cpbar=cp(2cp(x)-x)
% while imposing the interior dirichlet boundary conditions using the 
% direct approach, that is,
% u_left(ghost_left)=u_right(ghost_left) and
% u_right(ghost_right)=u_left(ghost_right).


% Include the cp_matrices folder (edit as appropriate)
addpath('../cp_matrices');
% add functions for finding the closest points
addpath('../surfaces');

%% Number of Iterations
t=cputime;
tic;
Nsteps=20; % Number of iterations
T=100*eps; % Time delay
% Generate the vector for errors
errvals_inf=zeros(1, Nsteps);
maxrelerr=zeros(1,Nsteps);
maxu=zeros(1,Nsteps);

%% Plot the composite domain 
R=1;
r=0.85; %% huge difference for r=0.85 and r=0.86
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
dx=0.1/16; % grid size
x1d=(-4-5*dx:dx:4+5*dx)';
y1d=(-1.5-5*dx:dx:1+5*dx)';
% Generate meshgrid
[x,y]=meshgrid(x1d,y1d);
% Make into vectors
X=x(:); Y=y(:);
% Plot the grid points
plot(X, Y, '.', 'color', 0.75*[1 1 1]);
hold on;
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
 
 
%% Closest Point representations
% CP points with respect to the rings 1,2,3,4,5 (numbered from left to right)
% [cp_x1, cp_y1,dist1,bdy1]=cpRing(X,Y,r,R,[-3,0]);
% plot(cp_x1,cp_y1, 'bo')   
% [cp_x2, cp_y2,dist2,bdy2]=cpRing(X,Y,r,R,[-1.5,-0.5]);
% plot(cp_x2,cp_y2, 'yo')
% [cp_x3, cp_y3,dist3,bdy3]=cpRing(X,Y,r,R,[0,0]);
% plot(cp_x3,cp_y3, 'ko')
% [cp_x4, cp_y4,dist4,bdy4]=cpRing(X,Y,r,R,[1.5,-0.5]);
% plot(cp_x4,cp_y4, 'go')
% [cp_x5, cp_y5,dist5,bdy5]=cpRing(X,Y,r,R,[3,0]);
% plot(cp_x5,cp_y5, 'ro')

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
% [cpbar_x1, cpbar_y1, dist6]=cpRing(2*cp_x1-X, 2*cp_y1-Y,...
%     r,R,[-3,0]);
% plot(cpbar_x1,cpbar_y1,'bo')
% [cpbar_x2, cpbar_y2, dist7]=cpRing(2*cp_x2-X, 2*cp_y2-Y,...
%     r,R,[-1.5,-0.5]);
% plot(cpbar_x2,cpbar_y2,'yo')
% [cpbar_x3, cpbar_y3, dist8]=cpRing(2*cp_x3-X, 2*cp_y3-Y,...
%    r,R,[0,0]);
%  plot(cpbar_x3,cpbar_y3,'ko')
% [cpbar_x4, cpbar_y4, dist9]=cpRing(2*cp_x4-X, 2*cp_y4-Y,...
%     r,R,[1.5,-0.5]);
% plot(cpbar_x4,cpbar_y4,'go')
% [cpbar_x5, cpbar_y5, dist10]=cpRing(2*cp_x5-X, 2*cp_y5-Y,...
%     r,R,[3, 0]);
% plot(cpbar_x5,cpbar_y5,'ro')
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


%% Assign true solution 

% u_true= @(x,y) sin(2*x)
u_true = @(x,y) (sin(4*x) + sin(5*y))/2;
% u_true=@(x,y) 1+0*x;
u1_true =u_true(cp_x1,cp_y1);
u2_true=u_true(cp_x2,cp_y2);
u3_true =u_true(cp_x3,cp_y3);
u4_true=u_true(cp_x4,cp_y4);
u5_true =u_true(cp_x5,cp_y5);


 %% f values
syms x y;
f=-diff(u_true,x,2)-diff(u_true,y,2);
f=matlabFunction(f);
% f= @(x,y) 4*sin(4*x);
f_1=f(cp_x1,cp_y1);
f_2=f(cp_x2,cp_y2);
f_3=f(cp_x3,cp_y3);
f_4=f(cp_x4,cp_y4);
f_5=f(cp_x5,cp_y5);
% f_1=zeros(size(X));
% f_2=zeros(size(X));
% f_3=zeros(size(X));
% f_4=zeros(size(X));
% f_5=zeros(size(X));
% Outer boundary condtions
u_Gamma_1=u1_true;
u_Gamma_2=u2_true;
u_Gamma_3=u3_true;
u_Gamma_4=u4_true;
u_Gamma_5=u5_true;

% u_Gamma_1=ones(size(X));
% u_Gamma_2=ones(size(X));
% u_Gamma_3=ones(size(X));
% u_Gamma_4=ones(size(X));
% u_Gamma_5=ones(size(X));
% 
% 
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
% ghost1_1 gives the set of all ghost points with respect ring 1
ghost1_1=(sqrt((cpbar_x1-X).^2+(cpbar_y1-Y).^2)>...
      sqrt((X-cp_x1).^2+(Y-cp_y1).^2));
% ghost2_1 gives the set of points within its negibouring ring
ghost2_1=((X+1.5).^2+(Y+0.5).^2<=R^2);
ghost3_1=((X+1.5).^2+(Y+0.5).^2>=r^2);
ring2=(ghost3_1 & ghost2_1);
% ghost4_1 and ghost5_1 gives the set of points near ring 1
ghost4_1=((X+3).^2+Y.^2<=R^2+2*dx);
ghost5_1=((X+3).^2+Y.^2>=r^2-2*dx);
%ghost_Gamma_1 gives the set of points we need to impose outer boundary
%condtion
ghost_Gamma_1=(~(ghost2_1 & ghost3_1) & ghost1_1 & ghost4_1 &ghost5_1 );
%ghost_12 gives the set of points we need to impose interior boundary
%condtion, i.e the boundary of disc 1 intersecting disc 2
ghost_12=(ghost1_1 & ring2 & ghost4_1 & ghost5_1);
%  plot(X(ghost_12),Y(ghost_12),'ko');
 plot(X(ghost_Gamma_1),Y(ghost_Gamma_1),'ro');

% ghost1_2 gives the set of all ghost points with respect disc 2
 ghost1_2=(sqrt((cpbar_x2-X).^2+(cpbar_y2-Y).^2)>...
       sqrt((X-cp_x2).^2+(Y-cp_y2).^2));
% ring1 and ring3 gives the set of points within its negibouring disc
 ghost2_2=((X+3).^2+Y.^2<=R^2);
 ghost3_2=((X+3).^2+Y.^2>=r^2);
 ring1 = (ghost2_2 & ghost3_2);
ghost4_2=(X.^2+Y.^2<=R^2);
ghost5_2=(X.^2+Y.^2>=r^2);
ring3=(ghost4_2 & ghost5_2);
% ghost6_2 and ghost7_2 gives the set of points near ring 2
ghost6_2=((X+1.5).^2+(Y+0.5).^2<=R^2+2*dx);
ghost7_2=((X+1.5).^2+(Y+0.5).^2>=r^2-2*dx);

%ghost_Gamma_2 gives the set of points we need to impose outer boundary
%condtion
ghost_Gamma_2=(~ring1 &~ring3 &ghost1_2& ghost6_2 & ghost7_2);
% %ghost_21 gives the set of points we need to impose interior boundary
% %condtion, i.e the boundary of disc 2 intersecting disc 2
 ghost_21=(ghost1_2 & ring1 & ghost6_2 & ghost7_2);
%ghost_23 gives the set of points we need to impose interior boundary
%condtion, i.e the boundary of disc 2 intersecting disc 3
ghost_23=(ghost1_2 & ring3 &ghost6_2 & ghost7_2);
% plot(X(ghost_Gamma_2),Y(ghost_Gamma_2),'ro');
% plot(X(ghost_21),Y(ghost_21),'ro')
% plot(X(ghost_23),Y(ghost_23),'ro')


% ghost1_3 gives the set of all ghost points with respect ring 3
ghost1_3=(sqrt((cpbar_x3-X).^2+(cpbar_y3-Y).^2)>...
    sqrt((X-cp_x3).^2+(Y-cp_y3).^2));
% ring4 give the set of points within its negibouring ring
ghost2_3=((X-1.5).^2+(Y+0.5).^2<=R.^2);
ghost3_3=((X-1.5).^2+(Y+0.5).^2>=r.^2);
ring4=(ghost2_3 & ghost3_3);

% ghost4_3 and ghost4_5 gives the set of points near ring 3
ghost4_3=(X.^2+Y.^2<=R.^2+2*dx);
ghost5_3=(X.^2+Y.^2>=r.^2-2*dx);
%ghost_Gamma_3 gives the set of points we need to impose outer boundary
%condtion
ghost_Gamma_3=(~ring2 & ~ring4 &ghost1_3& ghost4_3 &ghost5_3);
%ghost_32 gives the set of points we need to impose interior boundary
%condtion, i.e the boundary of disc 3 intersecting ring 2
ghost_32=(ghost1_3 & ring2 & ghost4_3&ghost5_3);
%ghost_34 gives the set of points we need to impose interior boundary
%condtion, i.e the boundary of disc 2 intersecting ring 4
ghost_34=(ghost1_3& ring4 & ghost4_3 &ghost5_3);
% plot(X(ghost_Gamma_3),Y(ghost_Gamma_3),'ro');
% plot(X(ghost_32),Y(ghost_32),'ro')

% ghost1_4 gives the set of all ghost points with respect to ring 4
ghost1_4=(sqrt((cpbar_x4-X).^2+(cpbar_y4-Y).^2)>...
    sqrt((X-cp_x4).^2+(Y-cp_y4).^2));
% ring5 gives the set of points within ring 5
ghost2_4=((X-3).^2+Y.^2<=R^2);
ghost3_4=((X-3).^2+Y.^2>=r^2);
ring5=(ghost2_4 & ghost3_4);
% ghost4_3 and ghost5_4 gives the set of points near ring 4
ghost4_4=((X-1.5).^2+(Y+0.5).^2<=R^2+2*dx);
ghost5_4=((X-1.5).^2+(Y+0.5).^2>=r^2-2*dx);
 %ghost_Gamma_4 gives the set of points we need to impose outer boundary
 %condtion
 ghost_Gamma_4=(~ ring3 & ~ring5 &ghost1_4& ghost4_4 & ghost5_4);
 %ghost_43 gives the set of points we need to impose interior boundary
 %condtion, i.e the boundary of ring 4 intersecting ring 3
ghost_43=(ghost1_4 & ring3 & ghost4_4 & ghost5_4);
%%ghost_45 gives the set of points we need to impose interior boundary
%condtion, i.e the boundary of ring 4 intersecting ring 5
ghost_45=(ghost1_4& ring5 & ghost4_4 & ghost5_4);
% plot(X(ghost_Gamma_4),Y(ghost_Gamma_4),'Bo');
% plot(X(ghost_43),Y(ghost_43),'ro')
% plot(X(ghost_45),Y(ghost_45),'yo')
 
% ghost1_5 gives the set of all ghost points with respect ring 5
 ghost1_5=(sqrt((cpbar_x5-X).^2+(cpbar_y5-Y).^2)>...
      sqrt((X-cp_x5).^2+(Y-cp_y5).^2));
 % ghost2_5 and ghost3_5 gives the set of points near disc 5
ghost2_5=((X-3).^2+Y.^2<=R^2+2*dx);
ghost3_5=((X-3).^2+Y.^2>=r^2-2*dx);
%ghost_Gamma_5 gives the set of points we need to impose outer boundary
%condtion
ghost_Gamma_5=(~ring4 & ghost1_5 &ghost2_5& ghost3_5);
%ghost_54 gives the set of points we need to impose interior boundary
%condtion, i.e the boundary of ring 5 intersecting ring 4
ghost_54=(ring4 & ghost1_5 & ghost2_5& ghost3_5);
% plot(X(ghost_54),Y(ghost_54),'Bo');
% plot(X(ghost_Gamma_5),Y(ghost_Gamma_5),'yo'); 

%% Building matrices used to impose dirichlet boundary conditions
 disp('building matrices to impose dirichlet boundary conditions')

% Interpolation matrices 
E12_bar=interp2_matrix(x1d,y1d,cpbar_x1(ghost_12),cpbar_y1(ghost_12),p);
E21_bar=interp2_matrix(x1d,y1d,cpbar_x2(ghost_21),cpbar_y2(ghost_21),p);
E23_bar=interp2_matrix(x1d,y1d,cpbar_x2(ghost_23),cpbar_y2(ghost_23),p);
E32_bar=interp2_matrix(x1d,y1d,cpbar_x3(ghost_32),cpbar_y3(ghost_32),p);
E34_bar=interp2_matrix(x1d,y1d,cpbar_x3(ghost_34),cpbar_y3(ghost_34),p);
E43_bar=interp2_matrix(x1d,y1d,cpbar_x4(ghost_43),cpbar_y4(ghost_43),p);
E45_bar=interp2_matrix(x1d,y1d,cpbar_x4(ghost_45),cpbar_y4(ghost_45),p);
E54_bar=interp2_matrix(x1d,y1d,cpbar_x5(ghost_54),cpbar_y5(ghost_54),p);
E_Gamma_1_bar=interp2_matrix(x1d,y1d,cpbar_x1(ghost_Gamma_1),cpbar_y1(ghost_Gamma_1),p);
E_Gamma_2_bar=interp2_matrix(x1d,y1d,cpbar_x2(ghost_Gamma_2),cpbar_y2(ghost_Gamma_2),p);
E_Gamma_3_bar=interp2_matrix(x1d,y1d,cpbar_x3(ghost_Gamma_3),cpbar_y3(ghost_Gamma_3),p);
E_Gamma_4_bar=interp2_matrix(x1d,y1d,cpbar_x4(ghost_Gamma_4),cpbar_y4(ghost_Gamma_4),p);
E_Gamma_5_bar=interp2_matrix(x1d,y1d,cpbar_x5(ghost_Gamma_5),cpbar_y5(ghost_Gamma_5),p);
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

M_ghost_Gamma_1=(I(ghost_Gamma_1,:)+E_Gamma_1_bar)/2;
M_ghost_Gamma_2=(I(ghost_Gamma_2,:)+E_Gamma_2_bar)/2;
M_ghost_Gamma_3=(I(ghost_Gamma_3,:)+E_Gamma_3_bar)/2;
M_ghost_Gamma_4=(I(ghost_Gamma_4,:)+E_Gamma_4_bar)/2;
M_ghost_Gamma_5=(I(ghost_Gamma_5,:)+E_Gamma_5_bar)/2;


M_1(ghost_12,:)=M_ghost_12;
M_2(ghost_21,:)=M_ghost_21;
M_2(ghost_23,:)=M_ghost_23;
M_3(ghost_32,:)=M_ghost_32;
M_3(ghost_34,:)=M_ghost_34;
M_4(ghost_43,:)=M_ghost_43;
M_4(ghost_45,:)=M_ghost_45;
M_5(ghost_54,:)=M_ghost_54;


M_1(ghost_Gamma_1,:)=M_ghost_Gamma_1;
M_2(ghost_Gamma_2,:)=M_ghost_Gamma_2;
M_3(ghost_Gamma_3,:)=M_ghost_Gamma_3;
M_4(ghost_Gamma_4,:)=M_ghost_Gamma_4;
M_5(ghost_Gamma_5,:)=M_ghost_Gamma_5;

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


E_ghost_Gamma_1=interp2_matrix(x1d, y1d, cp_x1(ghost_Gamma_1),cp_y1(ghost_Gamma_1),1);
E_ghost_Gamma_2=interp2_matrix(x1d, y1d, cp_x2(ghost_Gamma_2),cp_y2(ghost_Gamma_2),1);
E_ghost_Gamma_3=interp2_matrix(x1d, y1d, cp_x3(ghost_Gamma_3),cp_y3(ghost_Gamma_3),1);
E_ghost_Gamma_4=interp2_matrix(x1d, y1d, cp_x4(ghost_Gamma_4),cp_y4(ghost_Gamma_4),1);
E_ghost_Gamma_5=interp2_matrix(x1d, y1d, cp_x5(ghost_Gamma_5),cp_y5(ghost_Gamma_5),1);

disp ('done');

%% Building the Plotting matrices
disp ('building plotting matrices');
% 
tol=dx/2;
% Grid points on disc 1, 2,3,4 and 5
index_1=(sqrt((cpbar_x1-X).^2+(cpbar_y1-Y).^2)-...
       sqrt((X-cp_x1).^2+(Y-cp_y1).^2)<tol);
xx1=X(index_1); yy1=Y(index_1);
% Interpolate from the embedded space to the grid points on the 1st disc 
Eplot1=interp2_matrix(x1d,y1d,xx1,yy1, p);

index_2=(sqrt((cpbar_x2-X).^2+(cpbar_y2-Y).^2)-...
       sqrt((X-cp_x2).^2+(Y-cp_y2).^2)<tol);
xx2=X(index_2); yy2=Y(index_2);
% Interpolate from the embedded space to the grid points on the 2nd disc 
Eplot2=interp2_matrix(x1d,y1d,xx2,yy2, p);


index_3=(sqrt((cpbar_x3-X).^2+(cpbar_y3-Y).^2)-...
       sqrt((X-cp_x3).^2+(Y-cp_y3).^2)<tol);
xx3=X(index_3); yy3=Y(index_3);
% Interpolate from the embedded space to the grid points on the 3rd disc 
Eplot3=interp2_matrix(x1d,y1d,xx3,yy3, p);

index_4=(sqrt((cpbar_x4-X).^2+(cpbar_y4-Y).^2)-...
       sqrt((X-cp_x4).^2+(Y-cp_y4).^2)<tol);
xx4=X(index_4); yy4=Y(index_4);
% Interpolate from the embedded space to the grid points on the 4th disc 
Eplot4=interp2_matrix(x1d,y1d,xx4,yy4, p);

index_5=(sqrt((cpbar_x5-X).^2+(cpbar_y5-Y).^2)-...
       sqrt((X-cp_x5).^2+(Y-cp_y5).^2)<tol);
xx5=X(index_5); yy5=Y(index_5);
% Interpolate from the embedded space to the grid points on the 5th disc 
Eplot5=interp2_matrix(x1d,y1d,xx5,yy5, p);


% Grid points on disc 1, 2, 3, 4, 5 but not intersecting its neigbours
index_11=(index_1 & ~index_2);
xxx1=X(index_11); yyy1=Y(index_11);
Eplot_1=interp2_matrix(x1d,y1d,xxx1,yyy1, p);

index_22=(index_2 & ~index_1 & ~ index_3);
xxx2=X(index_22); yyy2=Y(index_22);
Eplot_2=interp2_matrix(x1d,y1d,xxx2,yyy2, p);

index_33=(index_3 & ~index_2 & ~ index_4);
xxx3=X(index_33); yyy3=Y(index_33);
Eplot_3=interp2_matrix(x1d,y1d,xxx3,yyy3, p);

index_44=(index_4 & ~index_3 & ~ index_5);
xxx4=X(index_44); yyy4=Y(index_44);
Eplot_4=interp2_matrix(x1d,y1d,xxx4,yyy4, p);

index_55=(index_5 & ~index_4);
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
    
    f_1(ghost_Gamma_1)=E_ghost_Gamma_1*u_Gamma_1;
    f_3(ghost_Gamma_3)=E_ghost_Gamma_3*u_Gamma_3;
    f_5(ghost_Gamma_5)=E_ghost_Gamma_5*u_Gamma_5;
    
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
f_2(ghost_Gamma_2)=E_ghost_Gamma_2*u_Gamma_2;
f_4(ghost_Gamma_4)=E_ghost_Gamma_4*u_Gamma_4;
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
title("Plot of u^{50} for CP Multiplicative Schwarz Method")

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
         [Eplot1*u1_true; Eplot_2*u2_true;Eplot3*u3_true;Eplot_4*u4_true;Eplot5*u5_true])));
     
     maxu(step)=max([Eplot1*u1; Eplot_2*u2;Eplot3*u3;Eplot_4*u4;Eplot5*u5]);
     maxrelerr(step) = errvals_inf(step)./ maxu(step);
     
 end
 figure(4);clf(4);
nn=1:Nsteps;
semilogy(nn,errvals_inf);
xlabel('Number of Iterations','FontSize',18); % x-axis label
ylabel('$$\|\tilde{v}^n-u\|_{\infty}$$','interpreter','latex','FontSize',18); % y-axis label
title ('Plot of maximal error vs number of iterations','FontSize',20);
set(gca,'FontSize',18)
 figure(5);clf(5);
nn=1:Nsteps;
semilogy(nn,maxrelerr);
xlabel('Number of Iterations', 'FontSize',18); % x-axis label
ylabel('$$\frac{\|\tilde{v}^n-u\|_{\infty}}{\|u\|_{\infty}}$$','interpreter','latex','FontSize',18); % y-axis label
title ('Plot of relative error vs number of iterations','FontSize',20);
set(gca,'FontSize',18)
toc
e=cputime-t;

