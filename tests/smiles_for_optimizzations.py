import cc_tools.gaussian_tools as gt
from cc_tools.utils import molecule_from_smiles
from cc_tools.pre_optimize import xtb_optimize_molecule
from rdkit.Chem import AllChem
import os

current_dir = os.path.dirname(os.path.abspath(__file__))


mols = [
    molecule_from_smiles(r"Cl[Pd-2]1([P+](C2=CC=CC=C2)(C3=CC=CC=C3)C4=CC=C5C(C=CC=C5)=C4C6=C([P+](C7=CC=CC=C7)1C8=CC=CC=C8)C=CC9=C6C=CC=C9)Cl", 
                         name="BINAP", charge=0, multiplicity=1),
    molecule_from_smiles(r"Cl[Pd-2]([P+](C1=CC=CC=C1)(C2=CC=CC=C2)C3=CC=CC=C3)([P+](C4=CC=CC=C4)(C5=CC=CC=C5)C6=CC=CC=C6)Cl", 
                         name="PPh3", charge=0, multiplicity=1)
]

cmd = "# opt freq b3lyp/def2svpp"

g_opt = gt.GaussianInputFile(file_directory=os.path.join(current_dir, "opt"))

for mol in mols:
    # AllChem.UFFOptimizeMolecule(mol) # not suitable for organometalic complex
    mol = xtb_optimize_molecule(mol)
    g_opt.molecule_to_gjf(mol, cmd=cmd)