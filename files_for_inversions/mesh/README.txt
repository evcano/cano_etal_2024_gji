The files of complete and trucated mesh are already included. However, we detail the steps to get them:

1. Run <mesh_circle.py complete_domain.yaml> or <mesh_circle.py truncated_domain.yaml>
	- This will use GMESH and save the mesh in <.msh> format

2. Run <LibGmsh2Specfem_convert_Gmsh_to_Specfem2D_official.py _domain.yaml -t A -b A -l A -r A>
	- This will convert the mesh from <.msh> to the format SPECFEM reads

3. Run <xconvert_external_layers_of_a_given_mesh_to_CPML_layers>
	- Manually add CPML thickness
	- add CPML layer at the top of mesh
	- CPML thickness is 125000 for complete mesh and 58000 for truncated mesh (roughly 3 elements)
