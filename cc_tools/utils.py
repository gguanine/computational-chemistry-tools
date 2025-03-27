from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
import cclib
import periodictable

import os

def xyz_block_to_cartesian(XYZBlock:str) -> str:
    # delete the first two lines of XYZ block to give Cartesian.
    return "\n".join(XYZBlock.splitlines()[2:]) 

def get_cartesian(mol:Chem.Mol) -> str:
    return xyz_block_to_cartesian(Chem.MolToXYZBlock(mol))
        
def set_cc_props(mol:Chem.Mol, name:str=None, charge=None, multiplicity=None):
    '''
    Set necessary molecular properties for computation chemistry: name, charge, and multiplicity\n
    Default properties will be calculated from molecule if not provided.
    '''
    
    if name is None:
        print("Molecule name is empty. Setting the name as InChI Key")
        name = Chem.inchi.MolToInchiKey(mol)  
    mol.SetProp("_Name", name) # Using _Name instead of Name to improve RDKit compatibility
    
    if charge is None:
        charge = Chem.GetFormalCharge(mol)
    mol.SetProp("Charge", str(charge))
        
    if multiplicity is None:
        print("Multiplicity is empty. Guess from structure")
        multiplicity = Descriptors.NumRadicalElectrons(mol) + 1
    if multiplicity != 1:
        print("Radical appears to be presented.")
    mol.SetProp("Multiplicity", str(multiplicity)) 
    
def print_mol_data(mol:Chem.Mol):
    # print molecule properties
    print("name =", mol.GetProp("_Name"))
    print("charge =", mol.GetProp("Charge"))
    print("multiplicity =", mol.GetProp("Multiplicity"))
    print(get_cartesian(mol), "\n\n") 
        
def molecule_from_smiles(smiles:str, name=None, charge=None, multiplicity=None) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=6)

    set_cc_props(mol, name, charge, multiplicity)
    print_mol_data(mol)
    
    return mol

def molecule_from_log(filename:str, name_suffix:str="out", allow_imag:bool = False) -> Chem.Mol:
    # When allow_imag is False, imaginary frequency is not allowed

    name = os.path.splitext(os.path.basename(filename))[0] + "_" + name_suffix
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

    # Obtain XYZ block of the converged structure
    atomlabel = [periodictable.elements[i].symbol for i in data.atomnos]
    coords = ["   ".join(map("{:.7f}".format, row)) for row in data.atomcoords[-1]] # the last coords (converged)
    atomNcoords = ["   ".join(row) for row in zip(atomlabel, coords)]
    cartesian = "\n".join(atomNcoords)
    XYZBlock = str(data.natom) + "\n\n" + cartesian

    # Embed molecule
    mol = Chem.MolFromXYZBlock(XYZBlock)
    set_cc_props(mol, name, data.charge, data.mult)
    print_mol_data(mol)
    print("Successfully load the molecule.")

    return mol        
        




    
    