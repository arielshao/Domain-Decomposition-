%% Plot the composite domain
% Plot the left square 
xleft1=-2;xleft2=0;
yleft1=-1;yleft2=1;
xleft = [xleft1, xleft2, xleft2, xleft1, xleft1];
yleft = [yleft1, yleft1, yleft2, yleft2, yleft1];

figure(1);clf;
h1=plot(xleft, yleft, 'b-', 'LineWidth', 2);
h1=fill(xleft,yleft,[0.5843 0.8157 0.9882]);
text(-1,0, '$$\tilde{\Omega}_{1}$$','Interpreter', 'latex','FontSize',26);
hold on;
% Plot the right square
xright1=0;xright2=2;
yright1=-1;yright2=1;
xright = [xright1, xright2, xright2, xright1, xright1];
yright = [yright1, yright1, yright2, yright2, yright1];
h2=plot(xright, yright, 'r-', 'LineWidth', 2);
h2=fill(xright,yright,[0.9882 0.5843 0.8157 ]);
text(1,0, '$$\tilde{\Omega}_{2}$$','Interpreter', 'latex','FontSize',26);
hold on;
xlim([-2.5 2.5]);
ylim([-1.5 1.5]);
axis equal
plot([xleft2+0.2,xleft2+0.2], [yleft1, yleft2],'-.b','LineWidth',2)
plot([xright1-0.2,xright1-0.2], [yright1, yright2],'-.r','LineWidth',2)

x1 = [0.2,0.55];
y1 = [0.76, 0.76];
a = annotation('textarrow',x1,y1);

text(-1,1.5, '$$\Omega_{1}$$','Interpreter', 'latex','FontSize',26);
x2 = [0.8,0.49];
y2 = [0.27, 0.27];
a = annotation('textarrow',x2,y2);

text(1,-1.5, '$$\Omega_{2}$$','Interpreter', 'latex','FontSize',26);
