infile = './infiles/PvE_LG9_GPP_summer_MC.txt';
nReps = 100;  %replicates of monte carlo sampling per site
master_out = 'MC_summary_PvE_LG9_GPP_summer.txt';
master_vt_out = 'MC_summary_vT_PvE_LG9_GPP_summer.txt';
fmaster = fopen(master_out,'w');
fmaster_vt = fopen(master_vt_out,'w');
CM_file = 'conf_matrix_counts_RProd_20191001_mod.txt';

addpath(genpath('/Users/brianhopkinson/Documents/MATLAB/Library/matGeom/matGeom')); %libarary for working with 3D meshes
addpath(genpath('/Users/brianhopkinson/Documents/MATLAB/Computer_Vision/pdollar_toolbox'));      %pdollar toolbox - various supporting functions - here imLabel

class_perm = [10 15 4 6 5 2 7 8 9 11 1 3 12 13 14]; %permute class order to desired: 
% [Algae, Galax, Orb, Past, Ssid, Apalm, Gvent, Antillo, Sea Rods, Rubble, not seen, Acerv, Sand, Unclass, Other]

nClasses = 15;

fid = fopen(infile,'r');
C = textscan(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\n');
outFiles   = C{1};
meshFiles  = C{2};
annFiles   = C{3};
irradFiles = C{4};
smartsFiles= C{5};
PvEFiles   = C{6};
RATE_DESIRED = C{7};
nProc = size(outFiles,1);

% load confusion matrix data and process to generate transition probabilty matrix
cm_perm = [14 6 15 9 8 7 4 3 5 1 10 11 13 12 2];
A = readmatrix(CM_file);
A = A(:,2:end); %strip off NaNs column from class labels
A = A'; %transpose so rows focus on predicted outcomes
A = A(cm_perm,cm_perm); %permute rows and columns so they correspond to nViewNet labels
TP = A./repmat(sum(A,2),1,size(A,2));  %transition probability matrix. given a prediction by nViewNet, the corresponding class row gives the probability that the given predictions should be have been the class corresponding to each column
%NOTE: this transition probablity matrix assumes relative frequency of classes in the prediction dataset match those in
%the annotation dataset.

[TPs, Itp] = sort(TP,2);
TPcumsum  = cumsum(TPs,2);

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
                         255, 165,  0;...  % 14 pink algae- orange               
                         ];                                 

ClassColorsList = ClassColorsList./255;


first = 1;
for i = 1:nProc
    
    %load data;
    %mesh data
    [V, F] = readMesh_off(meshFiles{i}); 
    nF = size(F,1);  %number of faces
    Mesh.V = V;
    Mesh.F = F;
    FArea = trimeshSurfaceArea_perFace(Mesh.V,Mesh.F);
    
    %annotation data
    meshAnns = importdata(annFiles{i});
    meshIdx = meshAnns(:,1)+1; %ZERO INDEXED FROM PYTHON!!! MATLAB IS 1 INDEXED!!!
    classID = meshAnns(:,2) + 1; % first class will used to indicate the mesh was not seen
    Pclass = zeros(size(F,1),1);
    Pclass(meshIdx) = classID;
    Pclass_aug = Pclass + 1; %increment - now 1 = not seen, all others class ids incremented as well.
    nClasses = 15;
    
    %load light field data
    load(irradFiles{i}); %E_global - total irradiance on each face (umol photons/m2/s), E_direct, E_diffuse, Time (time of day in hours), ix - indicator of whether face has line of sight to sun (i.e. recieving direct light).
    dur = ones(size(Time)); %duration between time points - right now all one hour but this could be made more flexible
    n_times = size(Time,1);
    
    %load PvE relationships
    fid = fopen(PvEFiles{i},'r');
    C = textscan(fid,'%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f');
    Pmax_avg = C{3};  % umol O2/m2/s
    Ek_avg   = C{4};  % umol photons/m2/s
    R_avg    = C{5};  %umol O2/m2/s
    Pmax_SEM  = C{6};   %Pmax standard error of mean (SEM)
    Ek_SEM    = C{7};  % Ek SEM
    R_SEM     = C{8};  %umol O2/m2/s
    
    %set up data receptacles
    ProdTot = zeros(nReps,1);
    ProdvT = zeros(nReps,n_times);
    ProdClass = zeros(nReps,nClasses);
    tic
    for j = 1:nReps

        %perturb PvE data from mean based on gaussian distribution 
        s = rng; 
        Pmax = Pmax_avg + Pmax_SEM.*randn(size(Pmax_SEM));
        Ek   = Ek_avg   + Ek_SEM  .*randn(size(Ek_SEM));
        R    = R_avg    + R_SEM   .*randn(size(R_SEM));

        if (strcmp(RATE_DESIRED{i}, 'GPP') == 1)
            R = zeros(size(R)); % make R zero for all classes for gross photosynthesis;
        end

        %convert units to mmol/m2/hr
        Pmax = Pmax * (3600/1000);
        R    = R    * (3600/1000);
        
        %perturb mesh element class predictions using transition probability matrix
        mesh_draw = rand(size(Pclass_aug)); %random uniform draws
        MC_class = zeros(nF,1);
        for k = 1:nF
           new_idx = find(TPcumsum(Pclass_aug(k),:) < mesh_draw(k));
           new_idx = new_idx(end) + 1; % find gave us the last element less than draw -> next one was selected by draw
           MC_class(k) = Itp(Pclass_aug(k), new_idx);
        end
        
        %create P_Ek_struct for simplified mesh PvE calculation
        Pmax_aug = [0; Pmax]; %augment with a zero for unclassified faces
        Ek_aug   = [1; Ek];   %augment with a one for unclassified faces - zero will give divide by zero errors
        R_aug    = [0; R];

        PvE_struct = zeros(nF, 3);
        PvE_struct(:,1) = Pmax_aug(MC_class);
        PvE_struct(:,2) = Ek_aug(MC_class);
        PvE_struct(:,3) = R_aug(MC_class);


        [PvT_byClass, PvT_reef] = PvE_on_mesh_func(Mesh,MC_class,PvE_struct, E_global,Time,nClasses,outFiles{i});
        PvT_byClass = PvT_byClass .* repmat(dur',nClasses,1);
        PvT_reef = PvT_reef .* dur';
        
        ProdClass(j,:) = sum(PvT_byClass,2)';
        ProdTot(j) = sum(PvT_reef);
        ProdvT(j,:) = PvT_reef;
   
    end  %end monte carlo reps
    toc
    %write to file
    save(outFiles{i},'ProdClass','ProdTot','ProdvT');
    
    %process MC data
    ProdTot_mean = mean(ProdTot);
    ProdTot_sd   = std(ProdTot);
    ProdClassF = ProdClass./repmat(ProdTot,1,nClasses);
    ProdClassF_mean = mean(ProdClassF,1);
    ProdClassF_sd   = std(ProdClassF,1);
    ProdvT_mean = mean(ProdvT,1);
    ProdvT_sd   = std(ProdvT,1);
    
    fprintf(fmaster,'%s\t',outFiles{i});
    fprintf(fmaster,'%f\t%f\t',ProdTot_mean, ProdTot_sd);
    for j =1:nClasses
        fprintf(fmaster,'%f\t%f\t',ProdClassF_mean(class_perm(j)), ProdClassF_sd(class_perm(j)));
    end
    fprintf(fmaster,'\n');
    
    for j = 1:n_times
        fprintf(fmaster_vt,'%f\t%f\t',ProdvT_mean(j),ProdvT_sd(j));
    end
    fprintf(fmaster_vt,'\n');
    
    

end

fclose(fmaster);
fclose(fmaster_vt);
