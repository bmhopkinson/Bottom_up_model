infile = './PvE_batch_trial.txt';

addpath(genpath('/Users/brianhopkinson/Dropbox/MATLAB/Library/matGeom/matGeom')); %libarary for working with 3D meshes
addpath(genpath('/Users/brianhopkinson/Dropbox/MATLAB/Computer_Vision/pdollar_toolbox'));      %pdollar toolbox - various supporting functions - here imLabel

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

Reef_PvE_outfile = 'Reef_PvE_fits.txt';
fek = fopen(Reef_PvE_outfile,'w');
first = 1;
for i = 1:nProc
    
    %load data;
    [V, F] = readMesh_off(meshFiles{i}); 
    nF = size(F,1);  %number of faces
    Mesh.V = V;
    Mesh.F = F;
    FArea = trimeshSurfaceArea_perFace(Mesh.V,Mesh.F);
    
    meshAnns = importdata(annFiles{i});
    meshIdx = meshAnns(:,1)+1; %ZERO INDEXED FROM PYTHON!!! MATLAB IS 1 INDEXED!!!
    classID = meshAnns(:,2) + 1; % class zero indicates no class was assigned
    Pclass = zeros(size(F,1),1);
    Pclass(meshIdx) = classID;
    nClasses = 15;
    
    
    %load PvE relationships
    fid = fopen(PvEFiles{i},'r');
    C = textscan(fid,'%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f');
    Pmax = C{3};  % umol O2/m2/s
    Ek   = C{4};  % umol photons/m2/s
    R    = C{5};  %umol O2/m2/s

    if (strcmp(RATE_DESIRED{i}, 'GPP') == 1)
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
    
    %load light field data
    load(irradFiles{i}); %E_global - total irradiance on each face (umol photons/m2/s), E_direct, E_diffuse, Time (time of day in hours), ix - indicator of whether face has line of sight to sun (i.e. recieving direct light).

    
    %determine total area of each class
    [path, name, ext] = fileparts(outFiles{i});
    areaFile = strcat(name,'_areas',ext);
    AreaPerClass = zeros(nClasses,1);
    for j = 1:nClasses
        AreaPerClass(j) = sum(FArea(Pclass == (j-1)));
    end
    AreaPerClass = AreaPerClass(class_perm);
    %dlmwrite(areaFile,AreaPerClass);

    [PvT_mesh, PvT_byClass, PvT_reef] = PvE_on_mesh_func(Mesh,Pclass_aug,PvE_struct, E_global,Time,nClasses,outFiles{i});
    %% output data and make some plots

    %plot legend first time through
    if first == 1 
        figure;
        ClassColorsPlot = reshape(flipud(ClassColorsList), nClasses, 1, 3);
        imagesc(ClassColorsPlot);
        set(gca,'XTick',[]); set(gca,'YTick',[]);
        imLabel(fliplr(Classes),'left');
        first = 0;
    end

    %write to file
    fout = fopen(outFiles{i},'w');
    fprintf(fout,'time\t PP_by_Class (mmol/m2/hr)\t');
    fprintf(fout,'%s\n',RATE_DESIRED{i});
    PvTout = PvT_byClass';
    nT = size(Time,1);
    for j = 1:nT
        fprintf(fout,'%f\t',[double(Time(j)) PvTout(j,:)]);
        fprintf(fout,'\n');
    end
    fclose(fout);
    
    %spatially resolved daily production
    Pint_mesh = sum(PvT_mesh,2);
    figure;
    colormap('hot');
    trisurf(Mesh.F, Mesh.V(:,1), Mesh.V(:,2), Mesh.V(:,3),Pint_mesh);
    %colorbar;
    axis image;
    shading flat;
    view([0,90]);
    
    %taxonomically resolved vs Time 
    figure;
    h = area(double(Time), PvT_byClass');
    xlabel('time (hrs)');
    ylabel('PP (mmol/m^{2}/hr)');

    %change face colors
    for j = 1:nClasses
        h(j).FaceColor =  ClassColorsList(j,:);
    end


    %PvE community 
    %incoming light data
    fid = fopen(smartsFiles{i},'r');
    C = textscan(fid,'%d %f %f %f %f'); 
    Time_smarts = C{1};
    H = C{2}.*(pi/180);% zenith angle - angle up from the horizon (H - 'height' angle) - 90 = on the horizon, 0 = directly overhead
    A = C{3}.*(pi/180); %azimuth angle - angle along the horizon  -zero degrees = north, increases in a clockwise fashion (e.g. 90 = east, 180 = south, 270 = west).
    E_norm = C{4}; %direct PAR irradiance normal to direction of sun angle (umol photons/m2/s)
    E_diff_horz = C{5}; %diffuse PAR on horizontal surface (umol photons/m2/s)
    fclose(fid);
    
    E_direct_horz = E_norm .* cos(H);
    E_tot = E_direct_horz + E_diff_horz;
    
    [E, iE] = sort(E_tot);
    NP     = PvT_reef(iE)';
    %make initial guess for parameter values
    mid = round(length(E)/2);
    beta0(1) = E(mid); %Ek guess
    beta0(2) =  mean(NP(end-1:end)); %Pmax - average last two values at highest irradiances
    if (strcmp(RATE_DESIRED{i}, 'NPP') == 1)
   	    beta0(3) = NP(1);   % respiration rate - y-offset
    end
    
    NPfit = [];
    Efit = linspace(E(1),E(end),100);
    %fit parameters of saturing exponential model with nonlinear fitting procedure.
    if (strcmp(RATE_DESIRED{i}, 'GPP') == 1)
        [betaFit, r, J, VCOV, mse] = nlinfit(E, NP, @satexp, beta0);
        NPfit = satexp(betaFit,Efit);
    else
        [betaFit, r, J, VCOV, mse] = nlinfit(E, NP, @satexp_yoff, beta0);
        NPfit = satexp_yoff(betaFit,Efit);
    end
    
    seb = sqrt(diag(VCOV));
    if (strcmp(RATE_DESIRED{i}, 'GPP') == 1)
        fprintf(1,'Ek: %6.4E +- %6.4E\t Pmax: %6.4E +- %6.4E\t Offset: %6.4E +- %6.4E\n',betaFit(1),seb(1),betaFit(2),seb(2), 0, 0);
        fprintf(fek,'mesh: %s\tEk:\t %6.4E\t%6.4E\tPmax:\t%6.4E\t%6.4E\tOffset:\t%6.4E\t%6.4E\n',outFiles{i},betaFit(1),seb(1),betaFit(2),seb(2), 0, 0);
    else
        fprintf(1,'Ek: %6.4E +- %6.4E\t Pmax: %6.4E +- %6.4E\t Offset: %6.4E +- %6.4E\n',betaFit(1),seb(1),betaFit(2),seb(2), betaFit(3), seb(3));
        fprintf(fek,'mesh: %s\tEk:\t %6.4E\t%6.4E\tPmax:\t%6.4E\t%6.4E\tOffset:\t%6.4E\t%6.4E\n',outFiles{i},betaFit(1),seb(1),betaFit(2),seb(2), betaFit(3), seb(3));
    end
    
    figure;
    plot(E_tot, PvT_reef,'ok');
    hold on;
    plot(Efit, NPfit);
    
    %areas
    planar_area = mesh_planar_area(Mesh.V,Mesh.F,'ConvexHull'); %currently convex hull of flattened mesh is used - this is not perfect.
    total_area  = sum(FArea);
 

end
fclose(fek);

