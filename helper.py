import os

def download_pdb_from_rcsb(pdb_id, data_path):
    input_protein_pdb_path = os.path.join(data_path, pdb_id+'.pdb')
    if os.path.isfile(input_protein_pdb_path) == False:
        os.system('wget https://files.rcsb.org/download/{} -P {}'.format(pdb_id, data_path))
    return input_protein_pdb_path