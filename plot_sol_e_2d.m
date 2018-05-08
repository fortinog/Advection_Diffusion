clear
nx = 101;
ny = 101;

a1 = dir('sol01_*.txt');
u = load(a1(1).name);
r0 = reshape(u(:,3),nx,ny);
x = reshape(u(:,1),nx,ny);
y = reshape(u(:,2),nx,ny);

for i = 1:1:length(a1)
    r = load(a1(i).name);
    r = reshape(r(:,3),nx,ny);
%     contour(x,y,r,linspace(-0.5,0.5,40))
   contour(x,y,r,50)

    %view(90,0)
    colorbar
    set(gca,'fontsize',18)
    drawnow

end
print -depsc2 solver_plot