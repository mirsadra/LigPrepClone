# tests/test_main.py

import unittest
from rdkit import Chem # type: ignore
from ligprepclone import main

class TestLigPrepClone(unittest.TestCase):
    def test_prepare_ligand(self):
        # Test with a simple molecule
        mol = Chem.MolFromSmiles('CCO')
        prepared_mol = main.prepare_ligand(mol)
        self.assertIsNotNone(prepared_mol)
        # Ensure the molecule has 3D coordinates
        conf = prepared_mol.GetConformer()
        self.assertEqual(conf.Is3D(), True)

if __name__ == '__main__':
    unittest.main()

