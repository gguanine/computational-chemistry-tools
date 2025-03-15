from rdkit import Chem
from rdkit.Chem import Descriptors

def xyz_block_to_cartesian(XYZBlock:str) -> str:
    # delete the first two lines of XYZ block to give Cartesian.
    return "\n".join(XYZBlock.splitlines()[2:]) 

def get_cartesian(mol:Chem.Mol) -> str:
    return xyz_block_to_cartesian(Chem.MolToXYZBlock(mol))

        
def set_cc_props(mol:Chem.Mol, name:str=None, charge:int=None, multiplicity:int=None):
    if name is None:
        print("Molecule name is empty. Setting the name as InChI Key")
        name = Chem.inchi.MolToInchiKey(mol)  
    mol.SetProp("_Name", name)
    
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
    # print molecule information
    print("name =", mol.GetProp("_Name"))
    print("charge =", mol.GetProp("Charge"))
    print("multiplicity =", mol.GetProp("Multiplicity"))
    print(get_cartesian(mol)) 
        
if __name__ == "__main__":        
    mol = Chem.MolFromSmiles("CCCCCCCO")
    set_cc_props(mol)

    print(mol.GetProp("_Name"))


    
    