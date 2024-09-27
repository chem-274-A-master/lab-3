#!/bin/bash

# Retrieve SDF for carbitol from PubChem
wget "https://pubchem.ncbi.nlm.nih.gov/rest/pug/conformers/00001FD200000005/SDF?response_type=save&response_basename=Conformer3D_COMPOUND_CID_8146" -O input_files/carbitol.sdf

# Make a mol2 file with charges calculated by bccc and atom types assigned by gaff2
cd antechamber
antechamber -i ../input_files/carbitol.sdf -fi sdf -o ../input_files/carbitol.mol2 -fo mol2 -c bcc -at gaff2

# Use antechamber to make a pdb to use with packmol
antechamber -i ../input_files/carbitol.mol2 -fi mol2 -o ../input_files/carbitol.pdb -fo pdb -c bcc -at gaff2
cd ../

# Make sure the pdb file has the correct residue names
sed -i 's/\bMOL\b/CAR/g' input_files/carbitol.pdb

# Make sure the mol2 file has the correct residue names
sed -i 's/\bMOL\b/CAR/g' input_files/carbitol.mol2

# Use parmchk to generate a parameter file
parmchk2 -i input_files/carbitol.mol2 -f mol2 -o input_files/carbitol.frcmod -s gaff2

# Make a carbitol library with the mol2 to load later
tleap -f input_files/carbitol_lib.in

# Use packmol to create a box of carbitol and water molecules
packmol < make_box.inp

# Get rid of chain identifiers in the pdb file
sed -i 's/ B /   /g' mixture.pdb

# Use tleap to create a prmtop and inpcrd file
tleap -f input_files/leap.in
