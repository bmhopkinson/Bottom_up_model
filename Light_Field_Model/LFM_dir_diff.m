%simple Light Field Model for reef productivity work 
% considers direct and diffuse light, but neglects scattering in the water or off the reef surfaces
% calculates incident light on 3D triangular mesh faces
% uses aabb bounding volume heirarchy to accelerate rayToSun calculation

addpath(genpath('/Users/brianhopkinson/Dropbox/MATLAB/Library/matGeom/matGeom')); %libarary for working with 3D meshes
addpath(genpath('/Users/brianhopkinson/Dropbox/MATLAB/Library/aabb-tree-master_2019')); %libarary for axis aligned BVH

% options
DO_LFM = 1;       %do the light field  model
DO_DISPLAY = 1;   %display results
DO_SMOOTHING = 0; %optional smoothing over time
DO_FLIP_Z    = 1; %flip z dimension? (441_simple = 1; 441 = 0,  A_646 = 0, B_647 =0, C_648 = 0, D_649 = 0) none of the metricRot meshes should need z_flip

%simulated diffuse data
diffuse_angles = [0 0.523 1.046 1.571]; %angle of surface to horizon in radians
K_d = -0.1; %diffuse PAR attenuation coeff

%data
meshFile = './meshes/0441_simple_2_model.off';
direct_irrad_File ='./light_data/smarts_results_KL_July.txt'; % from SMARTS model https://www.nrel.gov/rredc/smarts/
diffuse_irrad_File = './light_data/SMARTS_diffuse_coral_summer.txt';
outFile =  'LFM_dir_diff_simple_2.mat';

[V, F] = readMesh_off(meshFile);  
if(DO_FLIP_Z == 1)
    V(:,3) = -1*V(:,3);
end

[nF, ~] = size(F);

%calculate some mesh properties needed subseqeuntly
normals = meshFaceNormals(V, F); %unit vector normal to each face - used in determining how spread out direct light is on a face
centroids = meshFaceCentroids(V,F);

%determine angle of faces relative to horizontal for use in calculating direct irradiance intensity on each face
normal_sin_to_horiz = normals(:,3); % sin of angle between normal vector and horiz is z value of normal
normal_ang_to_horiz = asin(normal_sin_to_horiz);
face_ang_to_horiz = (pi/2) - abs(normal_ang_to_horiz);

E = meshEdges(F);  %mesh edges
E_len = meshEdgeLength(V, E, F); %length of mesh edges
p.rayOffset = 0.01*mean(E_len);  % distance by which ray origin is offset from face centroids when projecting back to sun to test for line of sight

meshMin = min(V);
meshMax = max(V);
distChar = norm(meshMax - meshMin); %characteristic distance of the mesh
p.rayLen = 5*distChar;

%setup dummy planar light sensor - near top of reef surface
[z_sort,idx] = sort(centroids(:,3),'descend');
sensor_z = mean(z_sort(1:round(0.1*size(z_sort,1)))); %mean of uppermost 10% of mesh elements as height of reef canopy
sensor_ang = 0; %flat sensor

%sun angle and direct irradiance data
fid = fopen(direct_irrad_File,'r');
C = textscan(fid,'%d %f %f %f %f'); 
Time = C{1};
n_time = size(Time,1);
Z_ang = C{2}; %solar zenith angle  - 0 directly overhead
Z_ang = Z_ang.*(pi/180);  %convert from degrees to radians
A_ang = C{3}; %azimuth angle - angle along the horizon  -zero degrees = north, increases in a clockwise fashion (e.g. 90 = east, 180 = south, 270 = west).
A_ang = A_ang.*(pi/180);
E_norm = C{4}; %direct PAR irradiance normal to direction of sun angle (umol photons/m2/s)
E_diff_horz = C{5}; %diffuse PAR on horizontal surface (umol photons/m2/s)
fclose(fid);
%correct for underwater
nw_na = 1.33; %ratio of refactive index in water/air - Kirk 2011 p 48
Z_ang_w = asin(sin(Z_ang)./nw_na);   % direction of light undewater due to refraction. via Kirk 2011 p 48
E_ang_w = pi/2 - Z_ang_w;% solar elevation angle - zero = on the horizon, 90 = directly overhead

%diffuse irradiance data
n_elev = 4;
n_azi  = 4;
fid =  fopen(diffuse_irrad_File,'r');
C = textscan(fid,['%f %f ', repmat('%f ',[1,14])]);
e_dif = (pi/180)* C{1,1}; %elevation angles of tilted surfaces, in radians
e_dif = e_dif([1:n_azi:size(e_dif,1)]);
a_dif = (pi/180)*C{1,2}; % azimuth angles of tilted surfaces, in radians
a_dif = a_dif(1:n_azi);
diffuse_irrad = {};

for i = 1:n_time
    this_dif = reshape(C{1,i+2},n_azi, n_elev); 
    diffuse_irrad{i} = this_dif';
    
end

if (DO_SMOOTHING == 1)
    nSmooth = 60;  %desired number of time points
    tsmooth   = linspace(Time(1),Time(end),nSmooth);
    Z_ang = interp1(Time, Z_ang, tsmooth)';
    E_ang_w = interp1(Time, E_ang_w, tsmooth)';
    A_ang = interp1(Time, A_ang, tsmooth)';
    E_norm = interp1(Time, E_norm,tsmooth)';
    E_diff_horz =  interp1(Time, E_diff_horz,tsmooth)';
    Time = tsmooth;
end

reflect = (1/2)*(sin(Z_ang - Z_ang_w).^2)./(sin(Z_ang + Z_ang_w).^2) + (1/2)*(tan(Z_ang - Z_ang_w).^2)./(tan(Z_ang + Z_ang_w).^2); %Fresnel's Eqn Kirk 2011 2.16 p 46
transmit = 1 - reflect;  %fraction of direct light transmitted through air-water interface


nTimes = size(Time,1);
%vector pointing toward the sun as viewed underwater
D = 1; %want unit length vector to represent sun angle
z = D.*sin(E_ang_w);
L = D.*cos(E_ang_w);
x = L.*sin(A_ang);
y = L.*cos(A_ang);

%flip vector
x = -1*x; y = -1*y; z = -1*z;
sunAngle = [x y z]; % (unit) vector describing angle of sun wrt local mesh reference frame (x, y, z)
assert(abs(norm(sunAngle(1,:)) - 1) < 1E-2, 'vector is not a unit vector');  %double check its a unit vector 

fprintf(1,'timepoint: ');
if(DO_LFM ==1)
    
    ix = zeros(nF, nTimes); %visibility indicator matrix - is there a line of sight between sun and mesh face
    fFull = zeros(nF, nTimes); %fraction of full direct irradiance being recieved by mesh face 
    E_diffuse = zeros(nF,nTimes);
    sensor_horiz = zeros(1,nTimes); %horizontal planar sensor
    sensor_follow = zeros(1,nTimes); %planar sensor that follows the sun angle - meant to simulate 4pi irradiance sensor (integrating)
    
    atten_dif = exp(K_d .* (-1*centroids(:,3))); 
    parfor i = 1:nTimes
        fprintf(1,'%d\t', i);
        fFull(:,i) = abs(dot(normals, repmat(sunAngle(i,:), nF, 1), 2));  %cos of angle between normal and sun vector - ignore sign b/c direction of normal is arbitrary. will need to test for intersection between face centroid and other faces in direction of sun angle
        ix(:,i) = rayToSun_BVH(V, F, sunAngle(i,:), p);
        X = repmat(a_dif',n_elev,1);
        Y = repmat(e_dif, 1, n_azi);
        E_diffuse(:,i) = atten_dif.*interp2(X,Y,diffuse_irrad{i},zeros(size(face_ang_to_horiz)),face_ang_to_horiz); % diffuse irradiance hitting face - right now azimuth is a dummy variable (zeros) but set so that it could be implemented later
        
        %sensor calculations
        sensor_horiz_dif(i) = exp(K_d .* (-1*sensor_z)) .* interp2(X,Y,diffuse_irrad{i},0,sensor_ang);
        sensor_follow_dif(i) = max(diffuse_irrad{i},[],'all');  %assume 4pi sensor detects max diffuse irrad.
    end
    d_dir = repmat(centroids(:,3),1,nTimes) ./repmat(sin(E_ang_w'),nF,1); %distance direct light must travel through water to face, account for increased pathlength at low solar elevation angles
    atten_dir =  exp(K_d .* -1*d_dir);
   
    E_direct = repmat(transmit'.*E_norm',nF,1).*atten_dir.* ix.*fFull;
    E_global = E_direct + E_diffuse;
    
    %sensor calculations
    d_sens = repmat(sensor_z,1,nTimes)./sin(E_ang_w');
    atten_sens = exp(K_d .* -1*d_sens);
    sensor_horiz_dir = transmit' .* atten_sens .* sin(E_ang_w') .* E_norm';
    sensor_follow_dir = transmit' .* atten_sens .* E_norm';
    sensor_horiz  = sensor_horiz_dir  + sensor_horiz_dif;
    sensor_follow = sensor_follow_dir + sensor_follow_dif;
    

    save(outFile,'E_direct','E_diffuse','E_global','Time','ix');
    save('sensor_data.mat','sensor_horiz','sensor_follow');

else
    load(outFile);
end


if(DO_DISPLAY == 1)
    for i = 1:nTimes
        plotMesh(F, V, E_global(:,i), sunAngle(i,:))
        pause(1)
    end
end
