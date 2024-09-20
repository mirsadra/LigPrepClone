# ligprepclone/main.py

import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
import sys
import os

def read_molecules(input_file):
    """
    Reads molecules from an input file.
    Supports SMILES (.smi), SDF (.sdf), and MOL (.mol) formats.
    """
    file_extension = os.path.splitext(input_file)[1].lower()
    if file_extension == '.smi':
        suppl = Chem.SmilesMolSupplier(input_file, titleLine=False)
    elif file_extension in ['.sdf', '.mol']:
        suppl = Chem.SDMolSupplier(input_file)
    else:
        print(f"Unsupported file format: {file_extension}")
        sys.exit(1)
    molecules = [mol for mol in suppl if mol is not None]
    return molecules

def prepare_ligand(mol):
    """
    Prepares a single ligand:
    - Adds hydrogens
    - Generates 3D coordinates
    - Optimizes geometry
    """
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol) == -1:
        print("Embedding failed.")
        return None
    if AllChem.UFFOptimizeMolecule(mol) != 0:
        print("Optimization failed.")
        return None
    mol = Chem.RemoveHs(mol)
    return mol

def write_molecules(molecules, output_file):
    """
    Writes molecules to an output file in SDF format.
    """
    writer = Chem.SDWriter(output_file)
    for mol in molecules:
        if mol is not None:
            writer.write(mol)
    writer.close()

def main():
    parser = argparse.ArgumentParser(description='Ligand Preparation Tool')
    parser.add_argument('-i', '--input', required=True, help='Input file (SMILES, SDF, MOL)')
    parser.add_argument('-o', '--output', required=True, help='Output file (SDF)')
    args = parser.parse_args()

    print("Reading input molecules...")
    molecules = read_molecules(args.input)
    print(f"Number of molecules read: {len(molecules)}")

    prepared_molecules = []
    for idx, mol in enumerate(molecules):
        print(f"Preparing molecule {idx+1}/{len(molecules)}")
        prepared_mol = prepare_ligand(mol)
        if prepared_mol is not None:
            prepared_molecules.append(prepared_mol)
        else:
            print(f"Molecule {idx+1} could not be prepared.")

    print("Writing prepared molecules to output file...")
    write_molecules(prepared_molecules, args.output)
    print("Ligand preparation completed successfully.")

if __name__ == '__main__':
    main()
