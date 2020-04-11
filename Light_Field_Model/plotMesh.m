function [ output_args ] = plotMesh( F, V, FColors, sunAngle )
%plot mesh data in figure 2. 

figure(1);
hold off
trisurf(F, V(:,1), V(:,2), V(:,3),FColors);
hold on
view([30 30]);
axis image;
shading flat;

M = max(V);
xmax = M(1);
ymax = M(2);

S = min(V);
xmin = S(1);
ymin = S(2);

x = linspace(xmin, xmax, 10);
y = linspace(ymin, ymax, 10);

[Xg, Yg] = meshgrid(x, y);
Xf = Xg(:);
Yf = Yg(:);

nf = size(Xf, 1);

quiver3(Xf, Yf, -1*ones(nf,1), sunAngle(1)*ones(nf,1), sunAngle(2)*ones(nf,1), sunAngle(3)*ones(nf,1));
colorbar;


end

