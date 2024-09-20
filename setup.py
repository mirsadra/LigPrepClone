# setup.py

from setuptools import setup, find_packages

setup(
    name='ligprepclone',
    version='0.1.0',
    description='Ligand preparation tool similar to Schrödinger\'s LigPrep',
    author='Mirsadra Molaei',
    packages=find_packages(),
    install_requires=[
        'rdkit'
    ],
)

