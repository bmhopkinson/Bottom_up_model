# Reef Productivity Bottom-up Model
These scripts calculate taxon-specific productivity on a 3D coral reef reconstruction as described as: 

Owen, D. P., M. H. Long, W. K. Fitt, B.M. Hopkinson. 2021. Taxon-specific primary production rates on coral reefs in the Florida Keys. Limnology and Oceanography. 66: 625-638 doi:10.1002/lno.11627.

# PvE (Productivity vs. Irradiance) Model Inputs
- a 3D triangular mesh reconstruction with class labels per face (e.g. Hopkinson et al PloSONE 15(3) e0230671)
- productivity rates per taxa as a function of irradiance 
- light intensities per mesh face over the course of a day (see Light Field Model below)

# Outputs
- the PvE model exports time-resolved, taxon-specific productivities and spatially resolved production rates
- errors can be estimated via Monte-Carlo resampling of the per taxa productivity estimates and probabilistic perturbation of face class labels



# Light Field Model

as described in Owen et al. 2021: 
The model considers direct and diffuse PAR separately and then sums the intensities of these two components to calculate the total PAR irradiance incident on each mesh face. PAR (400 – 700 nm) is not spectrally resolved. The direct and diffuse irradiance at the ocean surface was determined using the Simplified Model of Atmospheric Radiative Transfer of Sunshine (SMARTS Windows version 2.9.5i1.3) (Gueymard 2005) using the subtropical reference atmosphere, maritime aerosol model, a regional and tilted surface albedo of water/calm ocean, and a spectral range of 400 nm to 700 nm, for the hours of 6 AM to 7 PM on July 5 (“Summer”) and 7 AM to 5 PM on December 15 (“Winter”) at Little Grecian Reef. For site LG9, the SMARTS “Summer” date was shifted to June 27 to align with eddy covariance productivity data collected at this location (see below). The model was run on surfaces at four evenly spaced tilt and azimuth angles to parameterize diffuse irradiance as a function of these angles.  
As direct light passes through the air-water interface its intensity is reduced by reflection, which is treated using Fresnel’s equation, and its angle is modified by refraction (Kirk 2011). The intensity of direct light is attenuated with distance traveled through the water based on a diffuse attenuation coefficient for downwelling irradiance (Kd) of 0.1 m-1 (based on Zepp et al. 2008 and Ong et al. 2018), accounting for light absorption and scattering by the water column above the reef benthos. When direct light intercepts a mesh element, the angle of incidence is used to determine the intensity and a ray is projected from mesh element back to the sun to ensure the line of sight is not blocked by another part of the mesh. The line of sight test is accelerated using a bounding volume hierarchy. 
Diffuse PAR light is assumed to pass without loss through the air-water interface, but as it travels through the water it is attenuated based on a diffuse attenuation coefficient as described for direct irradiance. Diffuse light intensity striking a mesh face has angular dependence and this dependence was accounted for using the orientation of the mesh face and interpolating between outputs of the SMARTS model at discrete tilt and azimuth angles.  Direct and diffuse PAR irradiance incident on each face are then summed to determine the total PAR photosynthetic photon flux (µmol photons m-2 s-1) on each mesh face over time, which is then used in the “bottom up” model.
