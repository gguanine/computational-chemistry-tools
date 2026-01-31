from .config import HPC_NPROCS, HPC_MEM, GAUSSIAN_CHECK
from . import utils

from rdkit import Chem
from rdkit.Chem import AllChem

import os
  

def molecule_from_gjf(filename:str, name_suffix:str="_2") -> Chem.Mol:
    # a parser to read molecule from gaussian input file
    with open(filename, 'r') as f:
        text = f.read().strip()
        
    sections = text.split("\n\n")
        
    section2_lines = sections[2].splitlines()
    charge, multiplicity = map(int, section2_lines[0].split())  # the first line of section 2 in charge and mult
    
    natom = len(section2_lines[1:])  # the rest of section 2 is cartesian
    cartesian = "\n".join(section2_lines[1:])
    XYZBlock = str(natom) + '\n\n' + cartesian

    # embed molecule
    mol = Chem.MolFromXYZBlock(XYZBlock)
    
    name = sections[1].strip() # molecule name is always in the second section of .gjf file
    if name == '':
        print("Molecule name is empty. Setting the name as InChI Key")
        name = Chem.inchi.MolToInchiKey(mol) 
    print("Now reading", name)
        
    utils.set_cc_props(mol, name, charge, multiplicity)
    utils.print_mol_data(mol)
    print("Successfully load the molecule.")           
    
    return mol

class GaussianInputFile():
    def __init__(self, file_directory:str=None, nprocs:int=HPC_NPROCS, mem:str=HPC_MEM, check:bool=GAUSSIAN_CHECK):
        self.nprocs = nprocs
        self.mem = mem
        self.check = check
        print(f"\n\n\nInitializing Gaussian inputfile with nprocs={nprocs}, mem={mem}, and check={check}")

        if file_directory is None:
            self.file_path = os.getcwd()
        else:
            self.file_path = os.path.join(os.getcwd(), file_directory)
            if not os.path.exists(self.file_path):
                os.makedirs(self.file_path)
            print("Saving output files in directory", self.file_path)
            
    
    def molecule_to_gjf(self, mol:Chem.Mol, cmd:str, cmd2:str=""):
        if mol is None:
            print("**********Fail to load this molecule**********")
            return

        name = mol.GetProp("_Name")
        charge = mol.GetProp("Charge")
        multiplicity = mol.GetProp("Multiplicity")
        cartesian = utils.get_cartesian(mol)

        output_file_name = os.path.join(self.file_path, name + ".gjf")
        print("\n\n" + "-"*100 + "\nNow writing Gaussian input file " + output_file_name)
        print("with command", cmd)
        with open(output_file_name, "w+") as f:
            if self.nprocs is not None:
                f.write(f"%nprocshared={self.nprocs}\n")
            if self.mem is not None:
                f.write(f"%mem={self.mem}\n")
            if self.check is True:
                f.write(f"%chk={name}.chk\n")
            f.write(cmd + "\n\n")   
            f.write(name + "\n\n")
            f.write(f"{charge} {multiplicity}\n")
            f.write(cartesian) 
            f.write("\n\n" + cmd2 + "\n\n")
        print("done!\n" + "-"*100 + "\n")

def dir_logs_to_gjf(input_dir:str, output_dir:str, cmd:str, cmd2:str="", name_suffix:str=""):
    # finding  Gaussian output .log files in input direactory, convert them into gjf using command cmd.
    
    log_files = []
    for root, dirs, files in os.walk(input_dir):
        for file in files:
            if file.endswith('.log'):
                full_file_dir = os.path.join(root, file)
                print("Found:", full_file_dir)
                log_files.append(full_file_dir)
                
    mols = [utils.molecule_from_log(f, name_suffix=name_suffix) for f in log_files]            
                
    g_en = GaussianInputFile(file_directory=output_dir)
    for mol in mols:
        g_en.molecule_to_gjf(mol, cmd=cmd, cmd2=cmd2)

def dir_gjfs_to_gjf(input_dir:str, output_dir:str, cmd:str, cmd2:str="", MM_opt=False):
    # finding  Gaussian input .gjf files in input direactory, convert them into gjf using command cmd.
    
    original_files = []
    for root, dirs, files in os.walk(input_dir):
        for file in files:
            if file.endswith('.gjf'):
                full_file_dir = os.path.join(root, file)
                print("Found:", full_file_dir)
                original_files.append(full_file_dir)
                
    mols = [molecule_from_gjf(f, MM_opt=MM_opt) for f in original_files]            
                
    g_en = GaussianInputFile(file_directory=output_dir)
    for mol in mols:
        g_en.molecule_to_gjf(mol, cmd=cmd, cmd2=cmd2)

if __name__ == "__main__":
    pass






