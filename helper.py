import os

def download_pdb_from_rcsb(pdb_id, data_path):
    input_protein_pdb_path = os.path.join(data_path, pdb_id+'.pdb')
    if os.path.isfile(input_protein_pdb_path) == False:
        os.system('wget https://files.rcsb.org/download/{} -P {}'.format(pdb_id, data_path))
    return input_protein_pdb_path

def generate_pdb_4tleap(input_pdb, output_pdb):
    '''
    Source: https://github.com/choderalab/kinase-benchmark/issues/3
    '''
    infile = open(input_pdb, 'r')
    lines = infile.readlines()
    infile.close()

    outfile = open(output_pdb, 'w')
    for line in lines:
        if line[0:6] == 'ATOM  ':
            if line[13] != 'H': # might have this column wrong
                outfile.write(line)
        else:
            outfile.write(line)
    outfile.close()

def generate_tleap(input_pdb, prepi_file, frcmod_file,  tleap_in, topology_amber_prmtop, coordinate_amber_inpcrd):
    com_file = open(tleap_in,'w')
    com_file.write('''
    source leaprc.protein.ff14SB #Source leaprc file for ff14SB protein force field
    source leaprc.gaff #Source leaprc file for gaff
    source leaprc.water.tip3p #Source leaprc file for TIP3P water model
    loadamberprep {} #Load the prepi file for the ligand
    loadamberparams {} #Load the additional frcmod file for ligand
    mol = loadpdb {} #Load PDB file for protein-ligand complex
    solvatebox mol TIP3PBOX 8 #Solvate the complex with a cubic water box
    addions mol Cl- 0 #Add Cl- ions to neutralize the system
    saveamberparm mol {} {} #Save AMBER topology and coordinate files
    quit #Quit tleap program
    '''.format(prepi_file, frcmod_file, input_pdb, topology_amber_prmtop, coordinate_amber_inpcrd))
    com_file.close()