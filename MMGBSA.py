'''
Orignal source: https://github.com/quantaosun/Ambertools-OpenMM-MD/blob/682d775f18b8c321de019a0d943290b89756b0b2/Ambertools-OpenMM-MD_GBSA.ipynb
'''

import os
import argparse

from helper import count_ligands, create_mmpbsa_in, download_pdb_from_rcsb, generate_complex_tleapin, generate_pdb_4tleap, generate_protein_tleapin, mkdir_if_missing, simulation_openMM, split_ligands

parser = argparse.ArgumentParser(description='MMGBSA')
parser.add_argument(
    '--work_dir',
    default='./work_dir',
    type=str, help='Work directory')
parser.add_argument(
    '--data_folder',
    default='sample_data',
    type=str, help='Data folder name')
parser.add_argument(
    '--output_folder',
    default='output',
    type=str, help='Data folder name')
parser.add_argument(
    '--protein_pdb_id',
    default=None,
    type=str, help='None: if you use your custome pdf file; ID when you need to download, ex: 7L10')
parser.add_argument(
    '--protein_pdb_file',
    default='',
    type=str, help='Your pdb file in data_folder')
parser.add_argument(
    '--ligand_pdb_file',
    default='',
    type=str, help='Your pdb file in data_folder')
parser.add_argument(
    '--ligand_sdf_file',
    default=None,
    type=str, help='Your sdf file in data_folder')

# Arguments for Simulation
parser.add_argument('--rigidWater', action='store_true')
parser.add_argument('--no-rigidWater', dest='rigidWater', action='store_false')
parser.set_defaults(rigidWater=True)

parser.add_argument(
    '--ewaldErrorTolerance',
    default=0.0005,
    type=float, help='')
parser.add_argument(
    '--constraintTolerance',
    default=0.000001,
    type=float, help='')
parser.add_argument(
    '--simulation_step',
    default=1000000,
    type=int, help='')
parser.add_argument(
    '--equilibrationSteps',
    default=1000,
    type=int, help='')

parser.add_argument(
    '--simulation_platform',
    default='CPU',
    choices=['CPU', 'CUDA'], help='CPU/CUDA')

parser.add_argument(
    '--dcdReporter_step',
    default=10000,
    type=int, help='')

parser.add_argument(
    '--dataReporter_step',
    default=1000,
    type=int, help='')

parser.add_argument(
    '--checkpointReporter_step',
    default=10000,
    type=int, help='')

# Arguments for running ambertools and MMPBSA
parser.add_argument(
    '--AMBERHOME',
    default='/usr/local',
    type=str, help='The folder path contains amber.sh')

parser.add_argument(
    '--igb', choices=[1, 2, 5, 7, 8],
    default=5, help='')
parser.add_argument(
    '--number_frames_analysis',
    default=1,
    type=int, help='')
parser.add_argument(
    '--strip_mask',
    default='strip_mask=:WAT:Na+:Cl-:Mg+:K+',
    type=str, help='')
parser.add_argument(
    '--salt_concentration',
    default=0.15,
    type=float, help='')
args = parser.parse_args()


def protein_process(output_path, input_protein_pdb_path):
    # Fix your protein with PDBfixer with OpenMM Setup
    fix_protein_pdb = os.path.join(output_path, input_protein_pdb_path.split(
        '/')[-1].split('.')[-1] + '_fixed.pdb')
    os.system(
        'pdbfixer {} --output {}'.format(input_protein_pdb_path, fix_protein_pdb))

    # Fix protein pdb
    protein_amber_fix_pdb = os.path.join(
        output_path, 'protein.fixed.amber.pdb')
    os.system('pdbfixer {} --output {}'.format(fix_protein_pdb,
              protein_amber_fix_pdb))
    protein_amber_fix_reduceH_pdb = os.path.join(
        output_path, 'protein.fixed.amber_reduceH.pdb')
    os.system('reduce {} > {}'.format(
        protein_amber_fix_pdb, protein_amber_fix_reduceH_pdb))
    protein_amber_fix_reduceH_refix_pdb = os.path.join(
        output_path, 'protein.fixed.amber_reduceH.refix.pdb')
    os.system('pdb4amber -i {} -o {}'.format(protein_amber_fix_reduceH_pdb,
              protein_amber_fix_reduceH_refix_pdb))

    protein_tleap_pdb = os.path.join(output_path, 'protein_tleap.pdb')
    generate_pdb_4tleap(protein_amber_fix_reduceH_refix_pdb, protein_tleap_pdb)
    protein_tleap_in = os.path.join(output_path, 'protein.tleap.in')
    protein_prmtop = os.path.join(output_path, 'protein.prmtop')
    protein_inpcrd = os.path.join(output_path, 'protein.inpcrd')
    generate_protein_tleapin(
        protein_tleap_in, protein_tleap_pdb, protein_prmtop, protein_inpcrd)
    protein_tleap_out = os.path.join(output_path, 'protein.tleap.out')
    os.system('tleap -s -f {} > {}'.format(protein_tleap_in, protein_tleap_out))
    os.system('tail -n 500 {}'.format(protein_tleap_out))

    protein_prmtop_inpcrd = os.path.join(
        output_path, 'protein_prmtop_inpcrd.pdb')
    os.system('ambpdb -p {} -c  {} > {}'.format(protein_prmtop,
              protein_inpcrd, protein_prmtop_inpcrd))
    return protein_prmtop_inpcrd


def ligand_process(output_path, input_ligand_pdb_path):
    # Add hydrogens to ligand pdb
    ligand_wH_pdb = os.path.join(output_path, 'ligand_wH.pdb')
    os.system(
        'obabel -h -ipdb {} -opdb -O {}'.format(input_ligand_pdb_path, ligand_wH_pdb))

    # Convert to amber pdb for ligand
    ligand_amber_pdb = os.path.join(output_path, 'ligand_wH_4amber.pdb')
    os.system('pdb4amber -i {} -o {}'.format(ligand_wH_pdb, ligand_amber_pdb))

    # Generate ligand mol2 file
    ligand_amber_mol2 = os.path.join(output_path, 'ligand_wH_4amber.mol2')
    # ISSUE: change directory of sqm.out: https://github.com/choderalab/ambermini/issues/45
    # https://github.com/Amber-MD/amber-md.github.io/issues/7
    os.system('antechamber -fi pdb -i {} -fo mol2 -o {} -c bcc -pf y'.format(
        ligand_amber_pdb, ligand_amber_mol2))
    sqm_out_file = os.path.join(output_path, 'sqm.out')
    os.system('tail {}'.format(sqm_out_file))
    print('Make sure you see "--------- Calculation Completed ----------", otherwise something may have been wrong!')

    # Generate ligand prepi file
    ligand_prepi = os.path.join(output_path, 'ligand_wH_4amber_frommol2.prepi')
    os.system(
        'antechamber -i {} -fi mol2 -o {} -fo prepi -pf y'.format(ligand_amber_mol2, ligand_prepi))

    # Generate ligand frcmod file
    ligand_frcmod = os.path.join(
        output_path, 'ligand_wH_4amer__frommol2prepi.frcmod')
    os.system('parmchk2 -f prepi -i {} -o {}'.format(ligand_prepi, ligand_frcmod))

    return ligand_amber_pdb


def calculating_a_pair_protein_ligand(output_path, protein_prmtop_inpcrd, input_ligand_pdb_path):
    ligand_amber_pdb = ligand_process(output_path, input_ligand_pdb_path)
    # Combine protein_top_crd with ligand_amber into complex
    complex_pdb = os.path.join(output_path, 'complex.pdb')
    os.system('cat {} {} > {}'.format(
        protein_prmtop_inpcrd, ligand_amber_pdb, complex_pdb))

    # Convert complex pdb into Amber pdb
    complex_amber_pdb = os.path.join(output_path, 'complex.amber.pdb')
    os.system('pdb4amber {} > {}'.format(complex_pdb, complex_amber_pdb))

    # Fix complex amber pdb
    complex_amber_fix_pdb = os.path.join(
        output_path, 'complex.fixed.amber.pdb')
    os.system('pdbfixer {} --output {}'.format(complex_amber_pdb,
              complex_amber_fix_pdb))
    complex_amber_fix_reduceH_pdb = os.path.join(
        output_path, 'complex.fixed.amber_reduceH.pdb')
    os.system('reduce {} > {}'.format(
        complex_amber_fix_pdb, complex_amber_fix_reduceH_pdb))
    complex_amber_fix_reduceH_refix_pdb = os.path.join(
        output_path, 'complex.fixed.amber_reduceH.refix.pdb')
    os.system('pdb4amber -i {} -o {}'.format(complex_amber_fix_reduceH_pdb,
              complex_amber_fix_reduceH_refix_pdb))

    complex_tleap_pdb = os.path.join(output_path, 'complex_tleap.pdb')
    generate_pdb_4tleap(complex_amber_fix_reduceH_refix_pdb, complex_tleap_pdb)

    complex_tleap_in = os.path.join(output_path, 'complex.tleap.in')
    complex_topology_amber_prmtop = os.path.join(output_path, 'complex.prmtop')
    complex_coordinate_amber_inpcrd = os.path.join(
        output_path, 'complex.inpcrd')
    generate_complex_tleapin(complex_tleap_pdb, ligand_prepi, ligand_frcmod,
                             complex_tleap_in, complex_topology_amber_prmtop, complex_coordinate_amber_inpcrd)
    complex_tleap_out = os.path.join(output_path, 'complex.tleap.out')
    os.system('tleap -s -f {} > {}'.format(complex_tleap_in, complex_tleap_out))
    os.system('tail -n 500 {}'.format(complex_tleap_out))

    # Simulation by openMM
    trajectory_dcd_file, log_file, checkpointReporter_file = simulation_openMM(output_path, complex_topology_amber_prmtop, complex_coordinate_amber_inpcrd,
                                                                           args.rigidWater, args.ewaldErrorTolerance, args.constraintTolerance,
                                                                           args.simulation_step, args.equilibrationSteps, args.platform_name,
                                                                           args.dcdReporter_step, args.dataReporter_step, args.checkpointReporter_step)

    os.environ['AMBERHOME'] = args.AMBERHOME
    os.system('echo $AMBERHOME')
    os.system('source "$AMBERHOME/amber.sh"')
    mmpbsa_infile = os.path.join(output_path, 'mmpbsa.in')
    mbondi = create_mmpbsa_in(
        mmpbsa_infile, args.igb, args.number_frames_analysis, args.salt_concentration, args.strip_mask)

    Output_name = 'FINAL_RESULTS_MMPBSA.dat'
    final_mmpbsa = os.path.join(output_path, Output_name)

    trajectory_mdcrd = os.path.join(output_path, "trajectory.mdcrd")
    os.system("cpptraj -p {} -y {} -x {}".format(complex_topology_amber_prmtop, trajectory_dcd_file, trajectory_mdcrd))

    os.system('ante-MMPBSA.py  -p {} -c com.prmtop -r rec.prmtop -l ligand.prmtop -s {} -n :LIG --radii {}'.format(
        complex_topology_amber_prmtop, args.strip_mask, mbondi))
    # MMPBSA = "MMPBSA.py -O -i mmpbsa.in -o " + str(final_mmpbsa) +  ".dat -sp " + str(pdb_ref) + " -cp com.prmtop -rp rec.prmtop -lp ligand.prmtop -y "  + str(trajectory_dcd)
    os.system("MMPBSA.py -O -i {} -o {} -sp {} -cp com.prmtop -rp rec.prmtop -lp ligand.prmtop -y ".format(
        mmpbsa_infile, final_mmpbsa, complex_topology_amber_prmtop, trajectory_mdcrd))


def main():
    data_path = os.path.join(args.work_dir, args.data_folder)
    output_path = os.path.join(args.work_dir, args.output_folder)
    mkdir_if_missing(output_path)
    input_protein_pdb_path = os.path.join(data_path, args.protein_pdb_file)
    input_ligand_pdb_path = os.path.join(data_path, args.ligand_pdb_file)
    input_ligand_sdf_path = None
    if args.ligand_sdf_file is not None:
        input_ligand_sdf_path = os.path.join(data_path, args.ligand_sdf_file)
    if args.protein_pdb_id is not None:
        input_protein_pdb_path = download_pdb_from_rcsb(
            args.protein_pdb_id, data_path)

    protein_prmtop_inpcrd = protein_process(
        output_path, input_protein_pdb_path)

    if input_ligand_sdf_path is not None:
        if count_ligands(input_ligand_sdf_path) > 1:
            ligand_pdb_file_list = split_ligands(input_ligand_sdf_path)
            for input_ligand_pdb_path in ligand_pdb_file_list:
                ligand_name_i = input_ligand_pdb_path.split(
                    '/')[-1].split('.')[0]
                output_path_i = os.path.join(output_path, ligand_name_i)
                mkdir_if_missing(output_path_i)
                calculating_a_pair_protein_ligand(
                    output_path_i, protein_prmtop_inpcrd, input_ligand_pdb_path)
        else:
            input_ligand_pdb_path = input_ligand_sdf_path.replace(
                '.sdf', '.pdb')
            os.system('obabel {} -O {}'.format(input_ligand_sdf_path,
                      input_ligand_pdb_path))
            calculating_a_pair_protein_ligand(
                output_path, protein_prmtop_inpcrd, input_ligand_pdb_path)
    else:
        calculating_a_pair_protein_ligand(
            output_path, protein_prmtop_inpcrd, input_ligand_pdb_path)


if __name__ == "__main__":
    main()
