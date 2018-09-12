%% Plot the composite domain
%Plot the square 

xleft1=-2;xleft2=0;
yleft1=-1;yleft2=1;
xleft = [xleft1, xleft2, xleft2, xleft1, xleft1];
yleft = [yleft1, yleft1, yleft2, yleft2, yleft1];
figure(1);clf;
h1=plot(xleft, yleft, 'b-', 'LineWidth', 2);
h1=fill(xleft,yleft,[0 0 1]);
hold on;
plot([xleft2;xleft2], [yleft1;yleft2],'r-', 'LineWidth', 2);
hold on;
xlim([-3 2]);
ylim([-2 2]);
omega1="$$\Omega_{1} $$";
% text(-1.5, 0, omega1, 'Interpreter','latex','Fontsize', 24);
hold on;
gamma1="$$\gamma_1 $$";
% text(0.25, 0, gamma1, 'Interpreter','latex','Fontsize', 24)


%Plot the circle
thetas=linspace(0,2*pi,100)';
r=ones(size(thetas));
[xright,yright]=pol2cart(thetas,r);
xright=xright(:); 
yright=yright(:);
h2=plot(xright, yright, 'r-', 'LineWidth', 2);
h2=fill(xright,yright,[1 0 0]);
xlim([-3 2]);
ylim([-2 2]);
hold on;
% plot(xright(xright<0), yright(xright<0),'b-','LineWidth',2);
omega2="$$\Omega_{2} $$";
% text(0, 0, omega2, 'Interpreter','latex','Fontsize', 24);
% gamma2="$$\gamma_2 $$";
% text(-1.25, 0, gamma2, 'Interpreter','latex','Fontsize', 24);

x3=[xright(xright<0); xleft2; xleft2];
y3=[yright(xright<0); yleft1; yleft2];
h3=plot(x3,y3,'k-', 'LineWidth',2);
h3=fill(x3,y3,[1 0 1])
% omega="$$\Omega $$"
% text(0, 0, omega, 'Interpreter','latex','Fontsize', 24)

omega12="$$\Omega_{1,2} $$"
text(-0.5, 0, omega12, 'Interpreter','latex','Fontsize', 24)
