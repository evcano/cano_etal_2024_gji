Below, we explain how to generate the target structure model used in our manuscript.

1. Run the simulation in <files_for_inversions/specfem2d_complete_mesh>
	- This will save the spectral element mesh and an homogeneous velocity model in binary files.

2. Run <get_h5_model_slice.py>
	- This will save a slice of the CSEM velocity model in a <.npy> file.
	- This code needs a chunk of the CSEM velocity model. Contact eduardo.valerocano@kaust.edu.sa to request it.

3. Run <get_model_interpolator.py>
	- This will generate a SCIPY interpolator object.

4. Run <impose_vel_mode.py parfile_velocity.yaml>
	- This will write the CSEM slice into binary files that SPECFEM can read

5. Run <get_gaussian_field.py parfile_gaussian_field.yaml>
	- This will generate random velocity perturbations.

6. Run <perturb_vel_model.py parfile_perturbations.yaml>
	- This will impose the random velocity perturbations in the SPECFEM binary files containing the velocity model.

7. Copy the velocity model binary files into <files_for_inversions/specfem2d_complete_mesh/DATA> and run <xsmooth>
	- This will smooth the velocity model to remove heterogeneities smaller than the mesh resolution.

8. Run <extend_vel_model_on_pml_elements.py parfile_extend_model.yaml>
	- This will extend the velocity model into CPML elements, which is required to avoid numerical errors.


To get the target velocity model in the truncated mesh (inversions I and III of our manuscript), follow this steps:

1. Run the simulation in <files_for_inversions/specfem2d_truncated_mesh>
	- This will save the spectral element mesh and an homogeneous velocity model in binary files.

2. Run <get_truncated_model.py>
	- This will copy the target velocity model from the complete mesh to the truncated mesh.

3. Run <extend_vel_model_on_pml_elements.py parfile_extend_model_truncated.yaml>
	- This will extend the velocity model into CPML elements, which is required to avoid numerical errors.
