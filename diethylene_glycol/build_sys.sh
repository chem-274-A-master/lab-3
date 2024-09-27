#!/bin/bash

# Retrieve SDF for diethylene_glycol from PubChem
wget "https://pubchem.ncbi.nlm.nih.gov/rest/pug/conformers/00001FD200000005/SDF?response_type=save&response_basename=Conformer3D_COMPOUND_CID_8146" -O diethylene_glycol.sdf

# Make a mol2 file with charges calculated by bccc and atom types assigned by gaff2
cd antechamber
antechamber -i ../diethylene_glycol.sdf -fi sdf -o ../diethylene_glycol.mol2 -fo mol2 -c bcc -at gaff2

# Use antechamber to make a pdb to use with packmol
antechamber -i ../diethylene_glycol.mol2 -fi mol2 -o ../diethylene_glycol.pdb -fo pdb -c bcc -at gaff2
cd ../

# Make sure the pdb file has the correct residue names
sed -i 's/\bMOL\b/DEG/g' diethylene_glycol.pdb

# Make sure the mol2 file has the correct residue names
sed -i 's/\bMOL\b/DEG/g' diethylene_glycol.mol2

# Use parmchk to generate a parameter file
parmchk2 -i diethylene_glycol.mol2 -f mol2 -o diethylene_glycol.frcmod -s gaff2

# Make a diethylene_glycol library with the mol2 to load later
cd input_files
tleap -f make_lib.in

# Use packmol to create a box of diethylene_glycol and water molecules
packmol < make_box.inp
cd ../

# Get rid of chain identifiers in the pdb file
sed -i 's/ B /   /g' mixture.pdb

# Use tleap to create a prmtop and inpcrd file
tleap -f input_files/leap.in

# Copy the prmtop and inpcrd files to the correct directory
cp *.prmtop ../
cp *.inpcrd ../
