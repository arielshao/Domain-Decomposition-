addpath('../cp_matrices');
addpath('../surfaces');
%% Generate grid points around the composite domain
dx=0.05/2;
x1d = (-1.5:dx:1.5)';
y1d = (-1.5:dx:1.5)';
[x,y]=meshgrid(x1d,y1d);
% make into vectors
x=x(:); y=y(:);
cpf = @(x, y)cpAnnulus(x, y, 0.8, 1)
f=@(x,y) 0*x+1;
g=@(x,y) 0*x;
u= poisson_2d_banding(f,g,dx,x1d,y1d,cpf);
%% Plot the solution
hold on;
thetas=linspace(0,2*pi,100)';
r=ones(size(thetas));
[xright,yright]=pol2cart(thetas,r);
xright=xright(:); 
yright=yright(:);
plot(xright, yright, 'r-', 'LineWidth', 2);
xlim([x1d(1)-dx*5, x1d(end)+dx*5]);
ylim([y1d(1)-dx*5, y1d(end)+dx*5]);
axis equal
hold on;