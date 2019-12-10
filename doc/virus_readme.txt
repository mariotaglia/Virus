README Virus

1. Generation of protein input files

Virus accepts the protein structure in the file "kap.txt". To generate this file from a pdb, use pdb_to_coarse.m. This program generates a .txt file that should be renamed to "kap.txt" and a xyz file that can be open with VMD.

2. Rotation of protein input files.

If necessary, use gen_rotations.m to produce N rotations of the protein. Rotations reorient the dipole vector of the protein to align it to N points homogeneously distributed on the surface of a sphere that contains the protein.

3. Definitions file

The file DEFINITIONS.txt contains the input parameters for the calculation.
Important keywords:

-	delta: discretization length in nm 
-	wall: set to 0 (no wall) or 1 (wall)
-	sigmaq: if wall = 1, sets the surface charge density in units of |e|/nm^2
-	vtkflag:  set to 1 to save vtk files (use with care, these are very large files).
-	infile: initial guess, in general use 0
-	dimx, dimy, dimz: system size in delta units
-	dielS : Relative dielectric constant of protein
-	kaptype: set to 2. inmediatly after, input protein position and rotation matrix (is it more convenient to rotate the protein using the gen_rotations.m program).
-	fdisfromfile: set to 1 to read and use dissociation fractions of amino acids from "in- fdissaa.???.dat". Use to fix the dissociation fraction (i.e. no charge regulation).
K0fromfile: set to 1 to read and use bulk K0 from "in-K0.???.dat". Use to accelerate calculations. K0 is the effective equilibrium constant for an aa in bulk solution that should be used to guarantee f = 0.5 for pH = pKa. This value is independent of the presence of the wall, but can change with salt concentration.


4. Output files
 
-	fdissaa.* and fdisbulk.* : contain dissociation fraction of the aminoacids in the protein and in the bulk
-	K0.* : output K0 files
-	F_tot, F_mix?, F_eq: Contributions to the total free energy
-	Gmean.dat: Total free-energy vs pH
-	real_coords.txt: positions of the aminoacids in real space after protein translation/rotation.
-	sumq.dat : net charge of protein
-	poten and qprotT: electrostatic potential and protein charge as a function of position (x,y,z in vtk, averaged on yz plane in dat).
-	version: GIT version

