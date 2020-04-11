%calculate photosynthetic rates over time as a function of irradiance on labeled mesh
function [PvT_mesh, PvT_byClass, PvT_reef, Time] = PvE_on_mesh_func(Mesh,Pclass_aug,PvE_struct, E_global,Time, nClasses,outfile)

    %% do some preliminary processing
    
    %calculate the area of each face (m^2)
    FArea = trimeshSurfaceArea_perFace(Mesh.V,Mesh.F);
    nT = size(Time,1);
    nF = size(Mesh.F,1);  %number of faces

    %% calculate photosynthetic rates 

    %calculate photosynthetic rate per face over time (mmol/face/hr)%
    PvT_mesh = zeros(nF, nT); %preallocate space in memory
    for i=1:nF   %loop over all the faces in the mesh
        for j = 1:nT  %on each mesh face, loop over all the time points
            PvT_mesh(i, j) = PvE(E_global(i, j),PvE_struct(i,:)) .* FArea(i);  %calculate photosynthetic rate per unit area and multiplyby area  of the face
      
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
    planar_area = mesh_planar_area(Mesh.V,Mesh.F,'ConvexHull'); %currently convex hull of flattened mesh is used - this is not perfect.
    PvT_byClass = PvT_byClass./planar_area; %mmol/m2/hr
    PvT_reef = sum(PvT_byClass,1) ; %spatial and taxonomically integrated production rates in mmol/m2/hr
    PvT_mesh = PvT_mesh./repmat(FArea,1,nT);
    return
    
    

  





