%calculate photosynthetic rates over time as a function of irradiance on labeled mesh

%load libraries
addpath(genpath('/Users/brianhopkinson/Documents/MATLAB/Library/matGeom/matGeom')); %libarary for working with 3D meshes
addpath(genpath('/Users/brianhopkinson/Documents/MATLAB/Computer_Vision/pdollar_toolbox'));      %pdollar toolbox - various supporting functions - here imLabel

%output file
outFile = './output/215_ML_v3_PvE_Grad_sectors_20190313.txt';

RATE_DESIRED = 'GPP' ; % GPP or NPP

%define input files
% meshFile = './meshes_metricRot/0623_20171013_mesh_metricRot.off';
% annFile = './mesh_labelings/0623_20171013_mesh_predictions.txt';
% irradFile = './light_data/LFM_dir_diff_0623_20171013.mat';

meshFile = './meshes_metricRot/215_ML_v3_mesh_metricRot.off';
annFile = './mesh_labelings/215_ML_v3_predictions_v3.txt';
irradFile = './light_data/LFM_dir_diff_215_ML_v3_diffuse_water_20190306.mat';
smartsFile = './light_data/smarts_results_KL_July.txt';
PvEFile = 'PvE_data_by_class_nViewNet_summer.txt';

%% segment parameters
center = [12.5 19.0, 0]; % gradient location in ML_215_v3 %center point for bounding planes
%center = [4.3 21.5, 0]; % EC location in ML_215_v3 %center point for bounding planes

%relative to this center point we want planes parallel to z-axis radiating out at regular angles
theta_s = 0;
theta_e = 2*pi;
theta_steps = 12;

% radius bounds
r_b = [0 5.63];  % gradient footprint
%r_b = [0 4.35];  % EC footprint

thetas = linspace(theta_s, theta_e,theta_steps+1);
disk_area = (pi*r_b(2).^2) - (pi*r_b(1).^2);
seg_planar_areas = repmat((disk_area./theta_steps),theta_steps,1);


%% load data and do some preliminary processing
%load metric mesh (meters)
[V, F] = readMesh_off(meshFile);  
nF = size(F,1);  %number of faces


%calculate the area of each face (m^2)
FArea = trimeshSurfaceArea_perFace(V,F);
MeshArea = sum(FArea,1);

%load mesh annotations 
meshAnns = importdata(annFile);
meshIdx = meshAnns(:,1)+1; %ZERO INDEXED FROM PYTHON!!! MATLAB IS 1 INDEXED!!!
classID = meshAnns(:,2) + 1; % class zero indicates no class was assigned
Pclass = zeros(size(F,1),1);
Pclass(meshIdx) = classID;
nClasses = 15;

% Class Map - nViewNet
% 0 unclassified -  black
% 1 Apalm - brown    
% 2 Acerv - light brown  
% 3 Orbicella - medium green  
% 4 Siderastrea siderea-   yellow                         
% 5 Porites astreoides - lime green 
% 6 Gorgonian - purple  
% 7 Antillogorgia - red                         
% 8 Sea Rods- pink
% 9 algae -  green                         
% 10 rubble - medium blue                         
% 11 sand - light blue 
% 12 unclassified - dark grey
% 13 other - light grey
% 14 pink algae- for now same green as "algae"                  
     

% Class Map - OLD (sort of - really intermediate)
% 0 unclassified 
% 1 Acerv 
% 2 algae
% 3 Antillogorgia
% 4 Apalm
% 5 Gorgonia ventalina
% 6 Orbicella 
% 7 Porites astreoides
% 8 Plexaurella
% 9 rubble
% 10 sand              
% 11 Siderastrea siderea

% area of each class

AreaPerClass = zeros(nClasses,1);
for i = 1:nClasses
    AreaPerClass(i) = sum(FArea(Pclass == (i-1)));
end
dlmwrite('AreaPerClass.txt',AreaPerClass);


%displayClassifiedMesh(V,F, Pclass);

%load light field data
load(irradFile); %E_global - total irradiance on each face (umol photons/m2/s), E_direct, E_diffuse, Time (time of day in hours), ix - indicator of whether face has line of sight to sun (i.e. recieving direct light).
nT = size(Time,1);


%load PvE relationships
fid = fopen(PvEFile,'r');
C = textscan(fid,'%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f');
Pmax = C{3};  % umol O2/m2/s
Ek   = C{4};  % umol photons/m2/s
R    = C{5};  %umol O2/m2/s

if (strcmp(RATE_DESIRED, 'GPP') == 1)
	R = zeros(size(R)); % make R zero for all classes for gross photosynthesis;
end

%convert units to mmol/m2/hr
Pmax = Pmax * (3600/1000);
R    = R    * (3600/1000);

%create P_Ek_struct for simplified mesh PvE calculation
Pmax_aug = [0; Pmax]; %augment with a zero for unclassified faces
Ek_aug   = [1; Ek];   %augment with a one for unclassified faces - zero will give divide by zero errors
R_aug    = [0; R];
Pclass_aug = Pclass + 1; %increment - now 1 = unclassified, all others class ids incremented as well.
PvE_struct = zeros(nF, 3);
PvE_struct(:,1) = Pmax_aug(Pclass_aug);
PvE_struct(:,2) = Ek_aug(Pclass_aug);
PvE_struct(:,3) = R_aug(Pclass_aug);

%% calculate photosynthetic rates 

%calculate photosynthetic rate per face over time (mmol/face/hr)%
%PvE = @(E, Params) Params(1).*(1-exp(-E./Params(2))) - Params(3);  %PvE function Pmax*(1-exp(-E/Ek)) - R

PvT_mesh = zeros(nF, nT); %preallocate space in memory

for i=1:nF   %loop over all the faces in the mesh
    for j = 1:nT  %on each mesh face, loop over all the time points
        PvT_mesh(i, j) = PvE(E_global(i, j),PvE_struct(i,:)) .* FArea(i);  %calculate photosynthetic rate per unit area and multiply by area  of the face
    end %end loop on time (nT)
end %end loop on faces (nF);

%determine photosynthetic rates by class by summing all the faces in each class
PvT_byClass = zeros(nClasses,nT); %preallocate space in memory
for i = 1:nClasses
    for j = 1:nT
        PvT_byClass(i,j) = sum(PvT_mesh((Pclass_aug == i),j));
    end
end

%normalize rates per unit planar area
planar_area = mesh_planar_area(V,F,'ConvexHull'); %currently convex hull of flattened mesh is used - this is not perfect.
PvT_byClass = PvT_byClass./planar_area; %mmol/m2/hr
PvT_reef = sum(PvT_byClass,1) ; %spatial and taxonomically integrated production rates in mmol/m2/hr


% photosynthetic rates by segment
PvT_seg_cell= {};
PvT_class_allsegs = zeros(nClasses,nT);
for i = 1:theta_steps
    [Vc,Fc, indFaces] = sector_mesh(V,F, [thetas(i) thetas(i+1)],r_b, center);
    PvT_seg = PvT_mesh(indFaces,:);
    Pclass_seg = Pclass_aug(indFaces);
    
    PvT_seg_class = zeros(nClasses,nT);
    for j = 1:nClasses
        PvT_seg_class(j,:) = sum(PvT_seg((Pclass_seg == j),:));
    end
    PvT_class_allsegs = PvT_class_allsegs + PvT_seg_class;
    PvT_seg_cell{i} = PvT_seg_class./seg_planar_areas(i);
end

PvT_class_allsegs = PvT_class_allsegs./disk_area;


%% output data and make some plots

%write to file
fout = fopen(outFile,'w');
fprintf(fout,'time\t PP_by_Class (mmol/m2/hr)\t');
fprintf(fout,'%s\n',RATE_DESIRED);
for i  = 1:theta_steps
    fprintf(fout,'segment %d\n',i);
    PvTout = PvT_seg_cell{i};
    PvTout = PvTout';
    for j = 1:nT
        fprintf(fout,'%f\t',[double(Time(j)) PvTout(j,:)]);
        fprintf(fout,'\n');
    end
    fprintf(fout,'\n');
end
fclose(fout);


%PvsTime by class
figure(1);
h = area(double(Time), PvT_class_allsegs');
xlabel('time (hrs)');
ylabel('PP (mmol/m^{2}/hr)');

Classes = {'not seen','Apalm','Acerv','Orb','Ssid','Past', 'Gorg','Antill','Plexa','Algae',  'Rubble','Sand','Unclassified','Other','Pink algae' };
ClassColorsList = [  0,   0,   0;...  % 0 not_seen -  black
                         127,  51,   0;...  % 1 Apalm - brown    
                         160,  70,   0;...  % 2 Acerv - light brown  
                         100, 255, 25;...   % 3 Orbicella - medium green  
                         255, 255, 0 ;...   % 4 Siderastrea siderea-   yellow                         
                         201, 249, 138;...  % 5 Porites astreoides - lime green 
                         180,   0, 200;...  % 6 Gorgonian - purple  
                         255,   0,   0;...  % 7 Antillogorgia - red                         
                         255, 170, 238;...  % 8 Sea Rods- pink
                         50, 150, 50;...    % 9 algae -  green                         
                         90, 200, 255;...   % 10 rubble - medium blue                         
                         0, 255 , 255;...   % 11 sand - light blue 
                         100,100, 100;...   % 12 unclassified - dark grey
                         200,200, 200;...   % 13 other - light grey
                         50, 150, 50;...  % 14 pink algae- for now same green as "algae"                  
                         ];                                 
                     
ClassColorsList = ClassColorsList./255;

%change face colors
for i = 1:nClasses
    h(i).FaceColor =  ClassColorsList(i,:);
end

%legend
figure(2);
ClassColorsPlot = reshape(flipud(ClassColorsList), nClasses, 1, 3);
imagesc(ClassColorsPlot);
set(gca,'XTick',[]); set(gca,'YTick',[]);
imLabel(fliplr(Classes),'left');

% photosynthesis by sector around noon
noon_idx = 7;
P_byclass_noon = zeros(theta_steps, nClasses);
for i = 1:theta_steps
    thisPvT =PvT_seg_cell{i};
    P_byclass_noon(i,:) = thisPvT(:,noon_idx)';
    
end
figure(3);
b = bar(P_byclass_noon,'stacked');
xlabel('sector');
ylabel('PP (mmol/m^{2}/hr)');
jend = size(P_byclass_noon,2);
for i = 1:size(P_byclass_noon,2) %each class is considered a seperate object spread across the separate stacked bars
    b(i).FaceColor = 'flat';
    b(i).CData = ClassColorsList(i,:); 
end


%PvE community 
%incoming light data
% fid = fopen(smartsFile,'r');
% C = textscan(fid,'%d %f %f %f %f'); 
% Time_smarts = C{1};
% H = C{2}.*(pi/180);% zenith angle - angle up from the horizon (H - 'height' angle) - 90 = on the horizon, 0 = directly overhead
% A = C{3}.*(pi/180); %azimuth angle - angle along the horizon  -zero degrees = north, increases in a clockwise fashion (e.g. 90 = east, 180 = south, 270 = west).
% E_norm = C{4}; %direct PAR irradiance normal to direction of sun angle (umol photons/m2/s)
% E_diff_horz = C{5}; %diffuse PAR on horizontal surface (umol photons/m2/s)
% fclose(fid);
% 
% 
% E_direct_horz = E_norm .* cos(H);
% 
% figure(3)
% plot(E_direct_horz, PvT_reef,'ok');





