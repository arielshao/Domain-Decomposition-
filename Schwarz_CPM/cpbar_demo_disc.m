% adjust as appropriate

addpath('../cp_matrices');
addpath('../surfaces');
% Plot the disc
thetas=linspace(0,2*pi,100)';
r=ones(size(thetas));
[x,y]=pol2cart(thetas,r);
x=x(:); 
y=y(:);
h=plot(x, y);
hold on
axis equal
xlim([-2,2])
h=fill(x,y,[ 0.5843 0.8157 0.9882]);
% surface= '$$\Omega $$';
surface= '$$\mathcal{S} $$';
text(-0.5,0,surface, 'Interpreter','latex','Fontsize', 24);
xg=1.5; yg=1;
scatter(xg,yg,'ro', 'filled');
x_g= '$$ x_{g} $$';
text(xg+0.05,yg,x_g,'Interpreter','latex','Fontsize', 14)

[cpxg,cpyg,dist]=cpDisc(xg,yg);
[cpbarxg,cpbaryg,dist]=cpDisc(2*cpxg-xg,2*cpyg-yg);
scatter(cpxg,cpyg,'ro', 'filled');
cpx_g= '$$ cp(x_{g}) $$';
text(cpxg+0.05,cpyg,cpx_g,'Interpreter','latex','Fontsize', 14)
plot([xg,cpxg], [yg,cpyg],'-.');
scatter(cpbarxg, cpbaryg,'ro', 'filled')
cpbarx_g= '$$ 2cp(x_{g})-x_{g}=\bar{cp}(x_g)$$';

inequality = '$$ cp(x_{g})\neq\bar{cp}(x_g)$$';
text(0.5,1.2,inequality,'Interpreter','latex','Fontsize', 14);

text(cpbarxg+0.05,cpbaryg,cpbarx_g,'Interpreter','latex','Fontsize', 14)
plot([cpxg,cpbarxg], [cpyg,cpbaryg],'-.');

xx=-0.5; yy=-0.5;
x_gg= '$$ x=cp(x)=\bar{cp}(x)$$';
scatter(xx,yy,'ko', 'filled');
text(xx+0.05,yy,x_gg,'Interpreter','latex','Fontsize', 14)

