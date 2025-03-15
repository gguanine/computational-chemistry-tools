from .config import HPC_NPROCS, HPC_MEM, GAUSSIAN_CHECK
from . import utils

import periodictable
from rdkit import Chem
from rdkit.Chem import AllChem
import cclib

import os
  

def molecule_from_smiles(smiles:str, name=None, charge=None, multiplicity=None) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=1)
    AllChem.MMFFOptimizeMolecule(mol)

    utils.set_cc_props(mol, name, charge, multiplicity)
    utils.print_mol_data(mol)
    
    return mol

def molecule_from_log(filename:str, name_ending:str="-out", allow_imag:bool = False) -> Chem.Mol:
    # When allow_imag is False, imaginary frequency is not allowed

    name = os.path.splitext(os.path.basename(filename))[0] + "_" + name_ending
    print(f"\n\nNow reading molecule {name} in {filename}")

    data = cclib.io.ccopen(filename).parse()

    # Check if the Gaussian terminate nomally
    if not data.metadata["success"]:
        print("Error termination! Gaussian did not end normally.")
        return None

    # Check imaginary frequency
    try:
        if data.vibfreqs[0] < 0: # freq calculation may not perform
            print("Imaginary frequency is presented")
            if not allow_imag:
                print("***********Imaginary frequency is not allowed for this computation (by default).***********")
                return None
    except:
        pass


    XYZBlock = str(data.natom) + "\n\n"
    atomlabel = [periodictable.elements[i].symbol for i in data.atomnos]
    coords = ["   ".join(map("{:.7f}".format, row)) for row in data.atomcoords[-1]]
    atomNcoords = ["   ".join(row) for row in zip(atomlabel, coords)]
    cartesian = "\n".join(atomNcoords)
    XYZBlock += cartesian

    mol = Chem.MolFromXYZBlock(XYZBlock)
    mol.SetProp("Name", name)
    mol.SetProp("Charge", str(data.charge))
    mol.SetProp("Multiplicity", str(data.mult))
    mol.SetProp("Cartesian", cartesian) 
    #print_mol_data(mol)
    print("Successfully load the molecule.")

    return mol

def molecule_from_gjf(filename:str, output_ending_with:str="-out", MM_opt=False) -> Chem.Mol:
    # a parser to read molecule from gaussian input file
    with open(filename, 'r') as f:
        text = f.read().strip()
        sections = text.split("\n\n")
        
        mol_name = sections[1].strip()
        print("Now reading", mol_name)
        
        section2_lines = sections[2].splitlines()
        charge_multiplicity = section2_lines[0].split()
        charge = int(charge_multiplicity[0])
        multiplicity = int(charge_multiplicity[1])
        
        natom = len(section2_lines[1:])
        cartesian = "\n".join(section2_lines[1:])
        XYZBlock = str(natom) + '\n\n' + cartesian

        mol = Chem.MolFromXYZBlock(XYZBlock)
        
        if MM_opt is True:
            AllChem.MMFFOptimizeMolecule(mol)
            
        mol.SetProp("Name", mol_name)
        mol.SetProp("Charge", str(charge))
        mol.SetProp("Multiplicity", str(multiplicity))
        mol.SetProp("Cartesian", cartesian) 
        
        #print_mol_data(mol)
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
            f.write(f"%nprocshared={self.nprocs}\n")
            f.write(f"%mem={self.mem}\n")
            if self.check is True:
                f.write(f"%chk={name}.chk\n")
            f.write(cmd + "\n\n")   
            f.write(name + "\n\n")
            f.write(f"{charge} {multiplicity}\n")
            f.write(cartesian) 
            f.write("\n" + cmd2 + "\n\n")
        print("done!\n" + "-"*100 + "\n")

def dir_logs_to_gjf(input_dir:str, output_dir:str, cmd:str, name_ending=None):
    # finding  Gaussian output .log files in input direactory, convert them into gjf using command cmd.
    if name_ending is None:
        name_ending = input_dir
    
    cwd = os.getcwd()
    log_files = []
    for root, dirs, files in os.walk(os.path.join(cwd, input_dir)):
        for file in files:
            if file.endswith('.log'):
                full_file_dir = os.path.join(root, file)
                print("Found:", full_file_dir)
                log_files.append(full_file_dir)
                
    mols = [molecule_from_log(f, name_ending=name_ending) for f in log_files]            
                
    g_en = GaussianInputFile(file_directory=output_dir)
    for mol in mols:
        g_en.molecule_to_gjf(mol, cmd=cmd)

def dir_gjfs_to_gjf(input_dir:str, output_dir:str, cmd:str, MM_opt=False):
    # finding  Gaussian input .gjf files in input direactory, convert them into gjf using command cmd.
    
    cwd = os.getcwd()
    original_files = []
    for root, dirs, files in os.walk(os.path.join(cwd, input_dir)):
        for file in files:
            if file.endswith('.gjf'):
                full_file_dir = os.path.join(root, file)
                print("Found:", full_file_dir)
                original_files.append(full_file_dir)
                
    mols = [molecule_from_gjf(f, MM_opt=MM_opt) for f in original_files]            
                
    g_en = GaussianInputFile(file_directory=output_dir)
    for mol in mols:
        g_en.molecule_to_gjf(mol, cmd=cmd)

if __name__ == "__main__":
    pass






