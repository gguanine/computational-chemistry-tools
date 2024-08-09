from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
import cclib
import os

GAUSSIAN_NPROCS = 12
GAUSSIAN_MEM = "16GB"
GAUSSIAN_CHECK = True

at_num_symbol = \
    {1: 'H', 2: 'He',
     3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
     11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar',
     19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu',
     30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr',
     37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag',
     48: 'Cd', 49: 'In', 50: 'Sn', 51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe',
     55: 'Cs', 56: 'Ba', 57: 'La', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au',
     80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn',
     58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er',
     69: 'Tm', 70: 'Yb', 71: 'Lu'}

symbol_at_num = \
    {'H': 1, 'He':2,
     'Li': 3, 'Be':4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
     'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18,
     'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29,
     'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36,
     'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47,
     'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54,
     'Cs': 55, 'Ba': 56, 'La': 57, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79,
     'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86,
     'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68,
     'Tm': 69, 'Yb': 70, 'Lu': 71}

def XYZBlockToCartesian(XYZBlock:str) -> str:
    return "\n".join(XYZBlock.splitlines()[2:]) # delete the first two lines of XYZ block to give Cartesian.

def MolForDft_Smiles(smiles:str, name="", charge=None, multiplicity=None) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=1)
    AllChem.MMFFOptimizeMolecule(mol)

    print ("\n\n", "*"*50, sep="")

    # Set  molecule name
    if name == "":
        print("Molecule name is empty. Setting the name as InChI Key")
        name = Chem.inchi.MolToInchiKey(mol)  
    print(f"Now reading molecule:\nname = {name}")
    mol.SetProp("Name", name) 

    # Set molecule charge
    if charge is None:
        charge = Chem.GetFormalCharge(mol)
    print(f"charge = {charge}")
    mol.SetProp("Charge", str(charge)) 

    # Set molecule multiplicity
    if multiplicity is None:
        multiplicity = Descriptors.NumRadicalElectrons(mol) + 1
    if multiplicity != 1:
        print("Radical appears to be presented.")
    print(f"multiplicity = {multiplicity}")
    mol.SetProp("Multiplicity", str(multiplicity)) 

    # Set Cartesian
    cartesian = XYZBlockToCartesian(Chem.MolToXYZBlock(mol))
    print(cartesian)
    mol.SetProp("Cartesian", cartesian) 

    print ("*"*50)
    return mol

def MolForDft_Log(filename:str, name_ending:str="-out", allow_imag:bool = False) -> Chem.Mol:
    # When allow_imag is False, imaginary frequency is not allowed

    name = os.path.splitext(os.path.basename(filename))[0] + "_" + name_ending

    print(f"\n\nNow reading molecule {name} in {filename}")

    # Check if the Gaussian terminate nomally
    with open(filename, "r") as f:
        lines = f.readlines()
        for line in reversed(lines):
            if "Error termination" in line:
                print("Error termination! Gaussian did not end normally.")
                return None
            elif "Normal termination" in line:
                break

    data = cclib.io.ccopen(filename).parse()

    try:
        if data.vibfreqs[0] < 0: # freq calculation may not perform
            print("Imaginary frequency is presented")
            if not allow_imag:
                print("***********Imaginary frequency is not allowed for this computation (by default).***********")
                return None
    except:
        pass


    XYZBlock = str(data.natom) + "\n\n"
    atomlabel = [at_num_symbol[i] for i in data.atomnos]
    coords = ["   ".join(map("{:.7f}".format, row)) for row in data.atomcoords[-1]]
    atomNcoords = ["   ".join(row) for row in zip(atomlabel, coords)]
    cartesian = "\n".join(atomNcoords)
    XYZBlock += cartesian

    print(XYZBlock)
    print("charge = ", data.charge)
    print("multiplicity = ", data.mult)

    mol = Chem.MolFromXYZBlock(XYZBlock)
    mol.SetProp("Name", name)
    mol.SetProp("Charge", str(data.charge))
    mol.SetProp("Multiplicity", str(data.mult))
    mol.SetProp("Cartesian", cartesian) 
    print("Successfully load the molecule.")

    return mol

class GaussianInputFile():
    def __init__(self, file_directory:str=None, nprocs:int=GAUSSIAN_NPROCS, mem:str=GAUSSIAN_MEM, check:bool=GAUSSIAN_CHECK):
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
            
    
    def MolToGjf(self, mol:Chem.Mol, cmd:str, cmd2:str=""):
        if mol is None:
            print("**********Fail to load this molecule**********")
            return

        name = mol.GetProp("Name")
        charge = mol.GetProp("Charge")
        multiplicity = mol.GetProp("Multiplicity")
        cartesian = mol.GetProp("Cartesian")

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

def DirLogsToGjf(input_dir:str, output_dir:str, cmd:str, name_ending=None):
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
                
    mols = [MolForDft_Log(f, name_ending=name_ending) for f in log_files]            
                
    g_en = GaussianInputFile(file_directory=output_dir)
    for mol in mols:
        g_en.MolToGjf(mol, cmd=cmd)

if __name__ == "__main__":
    pass






