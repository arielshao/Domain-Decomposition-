% This code aims to implement the Multiplicative Schwarz Method to 
% solve the Poisson equation -\Delta u=f on 3 mutually overlapping rings
% with homogeneous or non-homoegenous dirichlet boundary conditions.

% This code deals with the outer dirichlet boundary condtions using 
% the extension approach, that is,
% u(x_g)=2u(cp(x_g))-u(cpbar(x_g)) where cpbar=cp(2cp(x)-x)
% while imposing the interior dirichlet boundary conditions using the 
% direct approach, that is,
% u_left(ghost_left)=u_right(ghost_left) and
% u_right(ghost_right)=u_left(ghost_right).

% ONLY right side interior bc is applied in this case

% We vary the inner radius to show that the effect of dimension of the annulus
% on the error.



% Include the cp_matrices folder (edit as appropriate)
addpath('../cp_matrices');
% add functions for finding the closest points
addpath('../surfaces');

%% Number of Iterations
t=cputime;
tic;
Nsteps=5; % Number of iterations
T=100*eps; % Time delay
% Generate the vector for errors
errvals_inf=zeros(1, Nsteps);
err=[];
%% Plot the composite domain 
R=1;
for r=0.05:0.05:0.95; 
LW='LineWidth';
% Plot the centre disc
thetas=linspace(0,2*pi,100)';
r_out=R*ones(size(thetas));
r_in=r*ones(size(thetas));
[x_centre_out,y_centre_out]=pol2cart(thetas,r_out);
[x_centre_in,y_centre_in]=pol2cart(thetas,r_in);
x_centre_out=x_centre_out(:); y_centre_out=y_centre_out(:);
x_centre_in=x_centre_in(:); y_centre_in=y_centre_in(:);
figure(1);clf;
plot(x_centre_out,y_centre_out, 'k-', LW,2);
hold on;
plot(x_centre_in,y_centre_in, 'k-', LW,2);
hold on;

x_left_out=x_centre_out-0.6;
y_left_out=y_centre_out+0.9;
x_left_in=x_centre_in-0.6;
y_left_in=y_centre_in+0.9;
plot(x_left_out,y_left_out, 'k-', LW,2);
plot(x_left_in,y_left_in, 'k-', LW,2);
hold on;

x_right_out=x_centre_out+0.6;
y_right_out=y_centre_out+0.9;
x_right_in=x_centre_in+0.6;
y_right_in=y_centre_in+0.9;
plot(x_right_out,y_right_out, 'k-', LW,2);
plot(x_right_in,y_right_in, 'k-', LW,2);
hold on;


xlim([-2 2])
ylim([-1.5 2.5])


 %% Generate grid points around the composite domain 
dx=0.025; % grid size
x1d=(-1.6-5*dx:dx:1.6+5*dx)';
y1d=(-1-5*dx:dx:1.9+5*dx)';
% Generate meshgrid
[x,y]=meshgrid(x1d,y1d);
% Make into vectors
X=x(:); Y=y(:);
% Plot the grid points
plot(X, Y, '.', 'color', 0.75*[1 1 1]);
hold on;
%% Initial conditions
u0= @(x,y) 1-sin(22*x).*sin(10*y);

% This makes the initial guess u0 into a vector
 u0 = u0(X,Y); 
%% Closest Point representations
% CP points with respect to the rings 1 (bottom),2 (right) and 3 (left)

[cp_x1, cp_y1,dist1,bdy1]=cpAnnulus(X,Y,r,R,[0,0]);
plot(cp_x1,cp_y1, 'bo')   
[cp_x2, cp_y2,dist2,bdy2]=cpAnnulus(X,Y,r,R,[0.6,0.9]);
plot(cp_x2,cp_y2, 'yo')
[cp_x3, cp_y3,dist3,bdy3]=cpAnnulus(X,Y,r,R,[-0.6,0.9]);
plot(cp_x3,cp_y3, 'ko')
% 
% Find the modified closest point representation of each grid point
% "cpbar" [Macdonald, Brandman, Ruuth 2011]: cpbar(x):=cp(2*cp(x)-x);

[cpbar_x1, cpbar_y1, dist4]=cpAnnulus(2*cp_x1-X, 2*cp_y1-Y,...
    r,R,[0,0]);
%   plot(cpbar_x1,cpbar_y1,'bo')
[cpbar_x2, cpbar_y2, dist5]=cpAnnulus(2*cp_x2-X, 2*cp_y2-Y,...
    r,R,[0.6,0.9]);
%  plot(cpbar_x2,cpbar_y2,'yo')
[cpbar_x3, cpbar_y3, dist6]=cpAnnulus(2*cp_x3-X, 2*cp_y3-Y,...
   r,R,[-0.6,0.9]);
%   plot(cpbar_x3,cpbar_y3,'ko')

 
%% Manufactured solution
uexact = @(x,y) (exp(sin(4*x)) + x.*sin(5*y))/2;
uexact_1 =uexact(cp_x1,cp_y1);
uexact_2=uexact(cp_x2,cp_y2);
uexact_3 =uexact(cp_x3,cp_y3);
f = @(x, y) 25 * x .* sin (5 * y) / 2 -8 * (-sin (4 * x) + cos (4 * x) .^ 2) .* exp (sin (4 * x));
f_1=f(cp_x1,cp_y1);
f_2=f(cp_x2,cp_y2);
f_3=f(cp_x3,cp_y3);

u_Gamma_1=uexact_1;
u_Gamma_2=uexact_2;
u_Gamma_3=uexact_3;
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
E1_1=interp2_matrix(x1d,y1d, cp_x1, cp_y1, 1);
E1_2=interp2_matrix(x1d,y1d, cp_x2, cp_y2, 1);
E1_3=interp2_matrix(x1d,y1d, cp_x3, cp_y3, 1);
I= speye(size(L));
% Modified Laplacian operator.
% Details can be found in [Macdonald & Brandman & Ruuth 2011]
L_1= E1_1*L-2*dim/dx^2*(I-E_1);
L_2= E1_2*L-2*dim/dx^2*(I-E_2);
L_3= E1_3*L-2*dim/dx^2*(I-E_3);
M_1=-L_1;
M_2=-L_2;
M_3=-L_3;
disp('done');
  
%% Find the ghost points we need for the Dirichlet boundary conditions 
dime=2;
order=2;
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));
% ring1, 2 and 3 give the set of points in ring1, 2 and 3
ring1=(X.^2+Y.^2<=R^2 & X.^2+Y.^2>=r^2);
ring2=((X-0.6).^2+(Y-0.9).^2<=R^2 & (X-0.6).^2+(Y-0.9).^2>=r^2);
ring3=((X+0.6).^2+(Y-0.9).^2<=R^2 & (X+0.6).^2+(Y-0.9).^2>=r^2);
% band1, 2 and 3 give the band around ring 1, 2 and 3 respectively
band1=(X.^2+Y.^2<=(R+bw*dx)^2 & X.^2+Y.^2>=(r-bw*dx)^2);
band2=((X-0.6).^2+(Y-0.9).^2<=(R+bw*dx)^2 & (X-0.6).^2+(Y-0.9).^2>=(r-bw*dx)^2);
band3=((X+0.6).^2+(Y-0.9).^2<=(R+bw*dx)^2 & (X+0.6).^2+(Y-0.9).^2>=(r-bw*dx)^2);
% ghost_Gamma_1, 2, 3 give the set of points we need to impose outer boundary
% condtion
ghost_Gamma_1=(~(ring1) & ~ring2 & band1);
ghost_Gamma_2=(~(ring2) & ~ring3 & band2);
ghost_Gamma_3=(~(ring3) & ~ring1 &  band3);

% The following ghost points are used to apply the interior dirichlet
% boundary conditions
ghost_12=(~(ring1) & ring2 & band1);
ghost_23=(~(ring2) & ring3 & band2);
ghost_31=(~(ring3) & ring1 & band3);

% TODO: check the positions of these points
% plot(X(ghost_Gamma_2),Y(ghost_Gamma_2),'ro')
% plot(X(ghost_31), Y(ghost_31),'yo')
% plot(X(ghost_32), Y(ghost_32),'bo')


%% Building matrices used to impose dirichlet boundary conditions
 disp('building matrices to impose dirichlet boundary conditions')

% Interpolation matrices 
E12_bar=interp2_matrix(x1d,y1d,cpbar_x1(ghost_12),cpbar_y1(ghost_12),p);
E23_bar=interp2_matrix(x1d,y1d,cpbar_x2(ghost_23),cpbar_y2(ghost_23),p);
E31_bar=interp2_matrix(x1d,y1d,cpbar_x3(ghost_31),cpbar_y3(ghost_31),p);

E_Gamma_1_bar=interp2_matrix(x1d,y1d,cpbar_x1(ghost_Gamma_1),cpbar_y1(ghost_Gamma_1),p);
E_Gamma_2_bar=interp2_matrix(x1d,y1d,cpbar_x2(ghost_Gamma_2),cpbar_y2(ghost_Gamma_2),p);
E_Gamma_3_bar=interp2_matrix(x1d,y1d,cpbar_x3(ghost_Gamma_3),cpbar_y3(ghost_Gamma_3),p);

%  M_ghost_12=(I(ghost_12,:)+E12_bar)/2;
%  M_ghost_23=(I(ghost_23,:)+E23_bar)/2;
%  M_ghost_31=(I(ghost_31,:)+E31_bar)/2;
% 
M_ghost_12=I(ghost_12,:);
M_ghost_23=I(ghost_23,:);
M_ghost_31=I(ghost_31,:);
% % 
M_ghost_Gamma_1=(I(ghost_Gamma_1,:)+E_Gamma_1_bar)/2;
M_ghost_Gamma_2=(I(ghost_Gamma_2,:)+E_Gamma_2_bar)/2;
M_ghost_Gamma_3=(I(ghost_Gamma_3,:)+E_Gamma_3_bar)/2;
% M_ghost_Gamma_1=I(ghost_Gamma_1,:);
% M_ghost_Gamma_2=I(ghost_Gamma_2,:);
% M_ghost_Gamma_3=I(ghost_Gamma_3,:);


M_1(ghost_12,:)=M_ghost_12;
M_2(ghost_23,:)=M_ghost_23;
M_3(ghost_31,:)=M_ghost_31;



M_1(ghost_Gamma_1,:)=M_ghost_Gamma_1;
M_2(ghost_Gamma_2,:)=M_ghost_Gamma_2;
M_3(ghost_Gamma_3,:)=M_ghost_Gamma_3;


% Bilinear interpolation from the grid points to the closest points of the
% % ghost_points (with respect to each ring)
E_ghost_12=interp2_matrix(x1d, y1d, cp_x1(ghost_12),cp_y1(ghost_12),1);
E_ghost_23=interp2_matrix(x1d, y1d, cp_x2(ghost_23),cp_y2(ghost_23),1);
E_ghost_31=interp2_matrix(x1d, y1d, cp_x3(ghost_31),cp_y3(ghost_31),1);

 
E_ghost_Gamma_1=interp2_matrix(x1d, y1d, cp_x1(ghost_Gamma_1),cp_y1(ghost_Gamma_1),1);
E_ghost_Gamma_2=interp2_matrix(x1d, y1d, cp_x2(ghost_Gamma_2),cp_y2(ghost_Gamma_2),1);
E_ghost_Gamma_3=interp2_matrix(x1d, y1d, cp_x3(ghost_Gamma_3),cp_y3(ghost_Gamma_3),1);
 disp ('done');

%% Building the Plotting matrices
disp ('building plotting matrices');

xx1=X(ring1); yy1=Y(ring1);
% Interpolate from the embedded space to the grid points on the 1st ring 
Eplot1=interp2_matrix(x1d,y1d,xx1,yy1, p);

xx2=X(ring2); yy2=Y(ring2);
% Interpolate from the embedded space to the grid points on the 2nd ring 
Eplot2=interp2_matrix(x1d,y1d,xx2,yy2, p);

xx3=X(ring3); yy3=Y(ring3);
% Interpolate from the embedded space to the grid points on the 3rd ring 
Eplot3=interp2_matrix(x1d,y1d,xx3,yy3, p);


xx12=X(ring1 & ring2 & ~ring3);yy12=Y(ring1 & ring2 & ~ring3);

xx23=X(ring2 & ring3 & ~ring1);yy23=Y(ring2& ring3 & ~ring1);
xx31=X(ring3 & ring1 & ~ring2);yy31=Y(ring3& ring1 & ~ring2);
 Eplot_12=interp2_matrix(x1d,y1d,xx12,yy12, p);
Eplot_23=interp2_matrix(x1d,y1d,xx23,yy23, p);
Eplot_31=interp2_matrix(x1d,y1d,xx31,yy31, p);

xxx1=X(ring1 & ~ring2 & ~ring3); yyy1=Y(ring1 & ~ring2 & ~ring3);
 Eplot_1=interp2_matrix(x1d,y1d,xxx1,yyy1, p);

xxx2=X(ring2 & ~ring1 & ~ring3); yyy2=Y(ring2 & ~ring1 & ~ring3);
Eplot_2=interp2_matrix(x1d,y1d,xxx2,yyy2, p);

xxx3=X(ring3 & ~ring1 & ~ring2); yyy3=Y(ring3 & ~ring1 & ~ring2);
Eplot_3=interp2_matrix(x1d,y1d,xxx3,yyy3, p);

% TODO: check the positions of these points
%plot(xxx3,yyy3,'go')
 disp ('done');
 
%% Plot initial solution u0
figure(2);
u0_1=Eplot1*u0;
u0_2=Eplot_2*u0;
u0_23=Eplot_23*u0;
u0_3=Eplot_3*u0;

plot2d_compdomain([u0_1;u0_2;u0_23;u0_3],[xx1;xxx2;xx23;xxx3],[yy1;yyy2;yy23;yyy3],dx,dx,2);
xlim([-2 2])
ylim([-1.5 2.5])
h=colorbar;
hold on;
plot(x_centre_out,y_centre_out, 'k-', LW,2);
hold on;
plot(x_centre_in,y_centre_in, 'k-', LW,2);
hold on;
plot(x_right_out,y_right_out, 'k-', LW,2);
hold on;
plot(x_right_in,y_right_in, 'k-', LW,2);
hold on;
plot(x_left_out,y_left_out, 'k-', LW,2);
hold on;
plot(x_left_in,y_left_in, 'k-', LW,2);

 %% Iterations
 counter=0;
 figure(3);
for step =1:Nsteps
     if step == 1
        u2=u0;
        u3=u0;
    end
     
 %Impose the dirichlet boundary conditions for ring1
%     f_1(ghost_12)=E_ghost_12*u2;
   f_1(ghost_12)=u2(ghost_12);
    f_1(ghost_Gamma_1)=E_ghost_Gamma_1*u_Gamma_1;
% f_1(ghost_Gamma_1)=u_Gamma_1(ghost_Gamma_1);
%  % Solve the pde in ring 1 using backslash
 u1=M_1\f_1;

%% Plot the solution 
uplot1=Eplot1*u1;
% extend the solution to the rest of the domain
 uplot_2=Eplot_2*u2;
 uplot_3=Eplot_3*u3;
 uplot_23=Eplot_23*u2;
% Plot the solution on the entire domain
plot2d_compdomain2([uplot1;uplot_2;uplot_23;uplot_3],[xx1;xxx2;xx23;xxx3],[yy1;yyy2;yy23;yyy3],dx,dx,3);
xlim([-2 2]);
ylim([-1.5 2.5]);
h=colorbar;
hold on;
plot(x_centre_out,y_centre_out, 'k-', LW,2);
hold on;
plot(x_centre_in,y_centre_in, 'k-', LW,2);
hold on;
plot(x_right_out,y_right_out, 'k-', LW,2);
hold on;
plot(x_right_in,y_right_in, 'k-', LW,2);
hold on;
plot(x_left_out,y_left_out, 'k-', LW,2);
hold on;
plot(x_left_in,y_left_in, 'k-', LW,2);

pause(T);

%Impose the dirichlet boundary conditions for ring2

%    f_2(ghost_23)=E_ghost_23*u3;
  f_2(ghost_23)=u3(ghost_23);
    f_2(ghost_Gamma_2)=E_ghost_Gamma_2*u_Gamma_2;
% f_2(ghost_Gamma_2)=u_Gamma_2(ghost_Gamma_2);

%  % Solve the pde in ring 1 using backslash
 u2=M_2\f_2;

%% Plot the solution 
uplot2=Eplot2*u2;
% extend the solution to the rest of the domain
 uplot_3=Eplot_3*u3;
 uplot_1=Eplot_1*u1;
 uplot_31=Eplot_31*u3;
% Plot the solution on the entire domain
plot2d_compdomain2([uplot2;uplot_3;uplot_31;uplot_1],[xx2;xxx3;xx31;xxx1],[yy2;yyy3;yy31;yyy1],dx,dx,3);
xlim([-2 2]);
ylim([-1.5 2.5]);
h=colorbar;
hold on;
plot(x_centre_out,y_centre_out, 'k-', LW,2);
hold on;
plot(x_centre_in,y_centre_in, 'k-', LW,2);
hold on;
plot(x_right_out,y_right_out, 'k-', LW,2);
hold on;
plot(x_right_in,y_right_in, 'k-', LW,2);
hold on;
plot(x_left_out,y_left_out, 'k-', LW,2);
hold on;
plot(x_left_in,y_left_in, 'k-', LW,2);

pause(T);

%Impose the dirichlet boundary conditions for ring 3
%   f_3(ghost_31)=E_ghost_31*u1;
   
  f_3(ghost_31)=u1(ghost_31);
    f_3(ghost_Gamma_3)=E_ghost_Gamma_3*u_Gamma_3;
% f_3(ghost_Gamma_3)=u_Gamma_3(ghost_Gamma_3);

%  % Solve the pde in ring 3 using backslash
 u3=M_3\f_3;

%% Plot the solution 
uplot3=Eplot3*u3;
% extend the solution to the rest of the domain
 uplot_1=Eplot_1*u1;
 uplot_2=Eplot_2*u2;
 uplot_12=Eplot_12*u1;
% Plot the solution on the entire domain
plot2d_compdomain2([uplot3;uplot_1;uplot_12;uplot_2],[xx3;xxx1;xx12;xxx2],[yy3;yyy1;yy12;yyy2],dx,dx,3);
xlim([-2 2])
ylim([-1.5 2.5])
h=colorbar;
hold on;
plot(x_centre_out,y_centre_out, 'k-', LW,2);
hold on;
plot(x_centre_in,y_centre_in, 'k-', LW,2);
hold on;
plot(x_right_out,y_right_out, 'k-', LW,2);
hold on;
plot(x_right_in,y_right_in, 'k-', LW,2);
hold on;
plot(x_left_out,y_left_out, 'k-', LW,2);
hold on;
plot(x_left_in,y_left_in, 'k-', LW,2);
pause(T);

counter=counter+1
%% The error calculated with respect to infinity norm
 errvals_inf(step)=max(max(abs([Eplot_1*u1; Eplot_12*u1;Eplot_2*u2; Eplot3*u3]-...
         [Eplot_1*uexact_1; Eplot_12*uexact_1;Eplot_2*uexact_2;Eplot3*uexact_3])));

end
err=[err errvals_inf(end)];
end

% 
%  figure(4);clf(4);
% nn=1:Nsteps;
% semilogy(nn,errvals_inf);
% xlabel('Number of Iterations'); % x-axis label
% ylabel('errvals_{inf}'); % y-axis label
% title ('Plot of errvals_{inf} vs number of iterations');
% toc
%%
figure(5);clf(5)
rr=0.05:0.05:0.95;
semilogy(rr,err);
xlabel('Inner disc radius r'); % x-axis label
ylabel('error after 5 iterations'); % y-axis label
title ('Plot of error vs r');
e=cputime-t;
% errvals_inf

