import os
import argparse

from helper import download_pdb_from_rcsb

parser = argparse.ArgumentParser(description='MMGBSA')
parser.add_argument(
        '--work_dir', 
        default='./work_dir',
        type=str, help='Work directory')
parser.add_argument(
        '--data_folder', 
        default='Sample_data',
        type=str, help='Data folder name')
parser.add_argument(
        '--protein_pdb_id', 
        default=None,
        type=str, help='None: if you use your custome pdf file; ID when you need to download, ex: 7L10')
parser.add_argument(
        '--protein_pdb_file', 
        default='',
        type=str, help='Your pdb file in data_folder (if exist)')
parser.add_argument(
        '--ligand_pdb_file', 
        default='',
        type=str, help='Your pdb file in data_folder (if exist)')
args = parser.parse_args()


def main():
    data_path = os.path.join(args.work_dir, args.data_folder)
    input_protein_pdb_path = os.path.join(data_path, args.protein_pdb_file)
    input_ligand_pdb_path = os.path.join(data_path, args.ligand_pdb_file)
    if args.protein_pdb_id is not None:
        input_protein_pdb_path = download_pdb_from_rcsb(args.protein_pdb_id, data_path)
    
    # Fix your protein with PDBfixer with OpenMM Setup
    fix_protein_pdb = os.path.join(data_path, input_protein_pdb_path.split('/')[-1].split('.')[-1] + '_fixed.pdb')
    os.system('pdbfixer {} --output {}'.format(input_protein_pdb_path, fix_protein_pdb))

    # Generate topolpogy and coordinate files by using gromacs
    coordinate_file = os.path.join(data_path, 'protein_processed.gro')
    topology_file = os.path.join(data_path,'topol.top')
    os.system('gmx pdb2gmx -f {} -o {} -ignh'.format(fix_protein_pdb, coordinate_file))

    # Combine topoly and coordinate into a new pdb file
    protein_top_crd_pdb = os.path.join(data_path, 'protein_top_crd.pdb')
    os.system('ambpdb -p {} -c {} > {}'.format(topology_file, coordinate_file, protein_top_crd_pdb))

    # Add hydrogens to ligand pdb
    ligand_wH_pdb = os.path.join(data_path, 'ligand_wH.pdb')
    os.system('obabel -h -ipdb {} -opdb -O {}'.format(input_ligand_pdb_path, ligand_wH_pdb))

    # Beautify pdb format for ligand
    # os.system('!gmx editconf -f {} -o $work_path/ligand2.gro'.format())

    # Convert to amber pdb for ligand
    ligand_amber_pdb = os.path.join(data_path, 'ligand_wH_4amber.pdb')
    os.system('pdb4amber -i {} -o {}'.format(ligand_wH_pdb, ligand_amber_pdb))

    # Generate ligand mol2 file
    ligand_amber_mol2 = os.path.join(data_path, 'ligand_wH_4amber.mol2')
    os.system('antechamber -fi pdb -i {} -fo mol2 -o {} -c bcc -pf y'.format(ligand_amber_pdb, ligand_amber_mol2))
    sqm_out_file = os.path.join(data_path, 'sqm.out')
    os.system('tail {}'.format(sqm_out_file))
    print('Make sure you see "--------- Calculation Completed ----------", otherwise something may have been wrong!')

    # Generate ligand prepi file
    ligand_prepi = os.path.join(data_path, 'ligand_wH_4amber_frommol2.prepi')
    os.system('antechamber -i {} -fi mol2 -o {} -fo prepi -pf y'.format(ligand_amber_mol2, ligand_prepi))

    # Generate ligand frcmod file
    ligand_frcmod = os.path.join(data_path, 'ligand_wH_4amer__frommol2prepi.frcmod')
    os.system('parmchk2 -f prepi -i {} -o {}'.format(ligand_prepi, ligand_frcmod))

    # Combine protein_top_crd with ligand_amber into complex
    complex_pdb = os.path.join(data_path, 'complex.pdb')
    os.system('cat {} {} > {}'.format(protein_top_crd_pdb, ligand_amber_pdb, complex_pdb))


if __name__ == "__main__":
    main()