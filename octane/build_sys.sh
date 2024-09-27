#!/bin/bash

# Retrieve SDF for octane from PubChem
wget "https://pubchem.ncbi.nlm.nih.gov/rest/pug/conformers/0000016400000001/SDF?response_type=save&response_basename=Conformer3D_COMPOUND_CID_356" -O octane.sdf

# Make a mol2 file with charges calculated by bccc and atom types assigned by gaff2
cd antechamber
antechamber -i ../octane.sdf -fi sdf -o ../octane.mol2 -fo mol2 -c bcc -at gaff2

# Use antechamber to make a pdb to use with packmol
antechamber -i ../octane.mol2 -fi mol2 -o ../octane.pdb -fo pdb -c bcc -at gaff2
cd ../

# Make sure the pdb file has the correct residue names
sed -i 's/\bMOL\b/OCT/g' octane.pdb

# Make sure the mol2 file has the correct residue names
sed -i 's/\bMOL\b/OCT/g' octane.mol2

# Use parmchk to generate a parameter file
parmchk2 -i octane.mol2 -f mol2 -o octane.frcmod -s gaff2

# Make a octane library with the mol2 to load later
cd input_files
tleap -f octane_lib.in

# Use packmol to create a box of octane and water molecules
packmol < make_box.inp
cd ../

# Get rid of chain identifiers in the pdb file
sed -i 's/ B /   /g' mixture.pdb

# Use tleap to create a prmtop and inpcrd file
tleap -f input_files/leap.in

# Copy the prmtop and inpcrd files to the correct directory
cp *.prmtop ../
cp *.inpcrd ../
