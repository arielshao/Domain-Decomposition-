function my_plot_new(XX,YY,ZZ,az,el,v, plot_range, col_range, counter, fix_axes)
surf(XX,YY,ZZ)
xlabel('x')
ylabel('y')
if fix_axes == 1
    axis([-1 1 -1 1 plot_range])
end
view(az-v*counter,el-0*counter)
shading(gca,'interp')
caxis(col_range)
end