function [Vc,Fc, indFaces] = sector_mesh(V,F,theta_b, r_b,center)
%slice a mesh (V,F) into a sector from the point 'center' between lower angle theta_l 
% and upper angle theta_ul; angles must be between 0 and 2pi (a bit confusingly 0 is the -x axis); 
% obtains all mesh faces within planes parallel to the z-axis radiating out at the defined angles
% returns the clipped mesh (Vc,Fc) and an index vector 


theta_l = theta_b(1);
theta_u = theta_b(2);

r_l = r_b(1);
r_u = r_b(2);

if (theta_l > 2*pi) || (theta_l < 0)
    error('theta_l must be between 0 and 2pi\n');
end

if (theta_u > 2*pi) || (theta_u < 0)
    error('theta_u must be between 0 and 2pi\n');
end

if (theta_l > theta_u)
    error('theta_l must be less than theta_u'); 
end


% account for center point
nV = size(V,1);
Vr = V - repmat(center,nV,1);

%convert to cylindrical coordinates 
r = ((Vr(:,1).^2 + Vr(:,2).^2)).^0.5;
theta = atan2(Vr(:,2), Vr(:,1)) + pi ;  %rescale angles between 0 and 2 pi (atan2 is -pi to +pi but that for a more confusing range).

%test if vertices are within upper and lower bounds of angles and radii
t_rl = (r > r_l);
t_ru = (r < r_u);
t_tl = (theta > theta_l);
t_tu = (theta < theta_u);
pass = (t_tl + t_tu + t_rl + t_ru) == 4;
indVertices = find(pass);


%extact vertices and faces that are within the bounds (pass the tests)
Vc = V(pass,:);

%remap vertex indices in the face structure - taken from geom3d library (clipMesh function)
refInds = zeros(size(indVertices));
for i = 1:length(indVertices)
    refInds(indVertices(i)) = i;
end

indFaces = sum(~ismember(F, indVertices), 2) == 0;
Fc = refInds(F(indFaces, :));

end

