%calculate photosynthetic rates over time as a function of irradiance on labeled mesh

%load libraries
addpath(genpath('/Users/brianhopkinson/Documents/MATLAB/Library/matGeom/matGeom')); %libarary for working with 3D meshes
addpath(genpath('/Users/brianhopkinson/Documents/MATLAB/Computer_Vision/pdollar_toolbox'));      %pdollar toolbox - various supporting functions - here imLabel

%output file
outFile = './output/215_ML_PvE_GPP_forcedLightField.txt';

RATE_DESIRED = 'GPP' ; % GPP or NPP
class_perm = [10 15 4 6 5 2 7 8 9 11 1 3 12 13 14]; %permute class order to desired: 
% [Algae, Galax, Orb, Past, Ssid, Apalm, Gvent, Antillo, Sea Rods, Rubble, not seen, Acerv, Sand, Unclass, Other]


%define input files
% meshFile = './meshes_metricRot/0623_20171013_mesh_metricRot.off';
% annFile = './mesh_labelings/0623_20171013_mesh_predictions.txt';
% irradFile = './light_data/LFM_dir_diff_0623_20171013.mat';

meshFile = './meshes_metricRot/215_ML_v3_mesh_metricRot.off';
annFile = './mesh_labelings/ML_215_v3_nViewNet_RProd_preds.txt';
PvEFile = './chamber_data/PvE_data_by_class_ReefProdPaper_20190929.txt';


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
E_forced = [88 450 911 1490 2059 2342 2609 2413 1991 1652 1343 529 202 23];
nT = size(E_forced,2);
Time = [6 7 8 9 10 11 12 13 14 15 16 17 18 19];

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
        PvT_mesh(i, j) = PvE(E_forced(j),PvE_struct(i,:)) .* FArea(i);  %calculate photosynthetic rate per unit area and multiply by area  of the face
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
PvT_byClass = PvT_byClass(class_perm,:); %permute class data (rows) to desired order
PvT_reef = sum(PvT_byClass,1) ; %spatial and taxonomically integrated production rates in mmol/m2/hr

%% output data and make some plots

%write to file
fout = fopen(outFile,'w');
fprintf(fout,'time\t PP_by_Class (mmol/m2/hr)\t');
fprintf(fout,'%s\n',RATE_DESIRED);
PvTout = PvT_byClass';
for i = 1:nT
    fprintf(fout,'%f\t',[double(Time(i)) PvTout(i,:)]);
    fprintf(fout,'\n');
end
fclose(fout);


%PvsTime by class

figure(1);
h = area(double(Time), fliplr(PvT_byClass'));
xlim([Time(1), Time(end)]);
ylim([-2 15]);
xlabel('time (hrs)');
ylabel('PP (mmol/m^{2}/hr)');

Classes = {'not seen','Apalm','Acerv','Orb','Ssid','Past', 'Gorg','Antill','Sea Rods','Algae',  'Rubble','Sand','Unclassified','Other','Galax' };
Classes = Classes(class_perm);
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
                         255, 165,  0;...  % 14 pink algae- orange               
                         ];                                 
                     
ClassColorsList = ClassColorsList./255;
ClassColorsList = ClassColorsList(class_perm,:);

%change face colors
ClassColorsList_inv = flipud(ClassColorsList);
for i = 1:nClasses
    h(i).FaceColor =  ClassColorsList_inv(i,:);
end

%legend
figure(2);
ClassColorsPlot = reshape(ClassColorsList, nClasses, 1, 3);
imagesc(ClassColorsPlot);
set(gca,'XTick',[]); set(gca,'YTick',[]);
imLabel(Classes,'left');

%PvE community 
figure(3)
plot(E_forced, PvT_reef,'ok');





