Below, we detail how to get the target noise source model.


1. Run the simulation in <files_for_inversions/specfem2d_complete_mesh>
	- This will save the spectral element mesh and an homogeneous velocity model in binary files.

2. Optional: Download maps of microseisms sources from https://sans.ethz.ch/

3. Run <get_sans_mean_model.py>
	- This will average maps downloaded from the SANS website
	- The average is stored in <mean_model.npy> and <mean_model_grd.npy>

4. Run <get_sans_model_interpolator.py>
	- This will generate a SCIPY interpolator object

5. Run <get_noise_distribution.py parfile_noise_model.yaml>
	- This will write the noise source distribution in SPECFEM binary files.

To get the noise source model in the truncated mesh (inversions I and III of our manuscript), follow this steps:

1. Run the simulation in <files_for_inversions/specfem2d_truncated_mesh>
	- This will save the spectral element mesh and an homogeneous velocity model in binary files.

2. Run <get_truncated_model.py>
	- This will copy the target noise source model from the complete mesh to the truncated mesh.

3. Run <mute_noise_model_on_pml_elements.py parfile_mute_model.yaml>
	- This will mute (set amplitude to zero) the noise source distribution in CPML boundary elements.
