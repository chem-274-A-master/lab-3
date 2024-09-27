#!/bin/bash

# Set variables for molecule name and abbreviation
MOLECULE_NAME="enflurane"
MOLECULE_ABBR="ENF"
COMPOUND_CID="3226"

# Retrieve SDF for the molecule from PubChem
wget "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${COMPOUND_CID}/record/SDF/?record_type=3d" -O ${MOLECULE_NAME}.sdf

# Make a mol2 file with charges calculated by bccc and atom types assigned by gaff2
cd antechamber
antechamber -i ../${MOLECULE_NAME}.sdf -fi sdf -o ../${MOLECULE_NAME}.mol2 -fo mol2 -c bcc -at gaff2

# Use antechamber to make a pdb to use with packmol
antechamber -i ../${MOLECULE_NAME}.mol2 -fi mol2 -o ../${MOLECULE_NAME}.pdb -fo pdb -c bcc -at gaff2
cd ../

# Make sure the pdb file has the correct residue names
sed -i "s/\bMOL\b/${MOLECULE_ABBR}/g" ${MOLECULE_NAME}.pdb

# Make sure the mol2 file has the correct residue names
sed -i "s/\bMOL\b/${MOLECULE_ABBR}/g" ${MOLECULE_NAME}.mol2

# Use parmchk to generate a parameter file
parmchk2 -i ${MOLECULE_NAME}.mol2 -f mol2 -o ${MOLECULE_NAME}.frcmod -s gaff2

# Make a molecule library with the mol2 to load later
cd input_files
sed "s/MOLECULE_NAME/${MOLECULE_NAME}/g; s/MOLECULE_ABBR/${MOLECULE_ABBR}/g" custom_lib.in.template > ../custom_lib.in

# Use packmol to create a box of molecule and water molecules
sed "s/MOLECULE_NAME/${MOLECULE_NAME}/g" make_box.inp.template > ../make_box.inp
cd ../

# Run tleap to create a custom library
tleap -f custom_lib.in

# Run packmol to create a box of molecules
packmol < make_box.inp

# Get rid of chain identifiers in the pdb file
sed -i 's/ B /   /g' mixture.pdb

# Use tleap to create a prmtop and inpcrd file
sed "s/MOLECULE_NAME/${MOLECULE_NAME}/g; s/MOLECULE_ABBR/${MOLECULE_ABBR}/g" input_files/leap.in.template > input_files/leap.in
tleap -f input_files/leap.in

# Copy the prmtop and inpcrd files to the correct directory
cp *.prmtop ../
cp *.inpcrd ../