# LigPrepClone

A ligand preparation tool that generates 3D coordinates and optimizes molecular geometries, similar to Schrödinger's LigPrep.

## Features

- Converts 2D molecular structures to 3D.
- Adds explicit hydrogens.
- Optimizes geometry using UFF.
- Supports SMILES, SDF, and MOL input formats.

## Installation

Ensure you have RDKit installed. It's recommended to use Anaconda:

```bash
conda create -c rdkit -n ligprepclone rdkit python=3.8
conda activate ligprepclone
