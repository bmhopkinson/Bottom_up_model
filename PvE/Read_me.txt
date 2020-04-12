These MATLAB functions are used to calculate photosynthesis as a function of irradiance over time.
Required datasets are: a 3D mesh, taxa labels for the 3D mesh,  photosynthesis as a function of irradiance for taxa, and light field data (e.g. produced by "Light_Field_Model"). 
1. Dependencies:  
	 - matGeom Toolbox (https://github.com/mattools/matGeom)
	 - pdollar Toolbox (https://pdollar.github.io/toolbox/)
	 install these toolboxes and modify paths in "PvE_on_mesh_batch.m" accordingly
	 
2. A trial dataset can be analyzed by running the PvE_on_mesh_batch.m script. 