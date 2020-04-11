function [ ix ] = rayToSun_BVH( V, F, sunAngle, p )
%test if ray from face to sun is blocked or clear
% V, F define triangluar mesh (V = vertices, F = faces);
% sunAngle is unit vector pointing from sun onto mesh (will need to invert to look back toward sun). 
% p = other paramters

%generate aabb BVH
bi = V(F(:,1),:); bj = V(F(:,1),:);
for ii = 2 : size(F,2)    
    bi = min(bi,V(F(:,ii),:)) ;
    bj = max(bj,V(F(:,ii),:)) ;
end

op.vtol = .67;
op.nobj = 32;

tic 
tr = maketree([bi,bj],op);
%fprintf(1,'done making tree\n');
%toc

%create rays to sun from faces
centroids = meshFaceCentroids(V,F);
raytosun = -1*sunAngle;  %invert sunAngle and then scale ray for BVH 
rayOrg = centroids + p.rayOffset*raytosun; %move ray origins slightly off centroid - otherwise all line of site rays will intersect with themselves
raytosun = p.rayLen * raytosun;
rays = [rayOrg repmat(raytosun, size(centroids,1), 1)];
tm = mapvert_ray(tr, rays);

%test for ray/mesh intersction using BVH acceleration
[ii, ir, it] = queryset(tr, tm, @rayInt_qs2, rays, F, V);
toc

nF = size(F, 1);

ix = ones(nF,1);  %indicates whether a face has a clear line to the sun, default is yes (1)
ix(ii) = 0;
    




end

