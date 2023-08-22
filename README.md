# Installation

Using with conda environment:
```
conda create -c conda-forge -c openbabel -c bioconda -c anaconda -n mmgbsa_env  python=3.7 openmm openmm-setup openbabel ambertools compilers pdbfixer gromacs babel rdkit
conda activate mmgbsa_env
```
Or:
```
```

# Usage
We can check the pipeline with sample data by using the command:
```
python MMGBSA.py --protein_pdb_file "protein.pdb" --ligand_pdb_file "ligand.pdb"
```

# OpenMM_MMGBSA

# FAQ
https://github.com/openmm/openmm/issues/2880
https://github.com/openmm/openmm/issues/2842
https://github.com/pablo-arantes/making-it-rain/issues/79
https://github.com/quantaosun/Ambertools-OpenMM-MD/issues/1

# Reference codebase
[Google colab](https://colab.research.google.com/drive/1OAF63N47PNpxuVR12RSqxEuznS6IjDM9#scrollTo=H-oBjHKEBPJY)