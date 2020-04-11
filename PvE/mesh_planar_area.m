function [planar_area] = mesh_planar_area(V,F,Method)

%calculate planar area of a mesh for normalization of productivity data.

if(Method == 'ConvexHull')
    Vflat = V(:,1:2);
    Fchull = convexHull(Vflat); %geom2d wrapper to matlab convhull
    planar_area = polygonArea(Fchull);
    %drawPolygon(Vflat, Fchull);  
else
    error('no such method to calculate mesh planar area');
end



end