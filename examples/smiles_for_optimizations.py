import cc_tools.gaussian_tools as gt
from cc_tools.utils import molecule_from_smiles
from cc_tools.pre_optimize import xtb_optimize_molecule
import os

current_dir = os.path.dirname(os.path.abspath(__file__))


mols = [
    molecule_from_smiles(r"[N+]1([Li-]2)=C(C3=[N+]2C=CC=C3)C=CC=C1", 
                         name="Libpy", charge=1, multiplicity=1)
]

cmd = "# opt freq hf/3-21g"

g_opt = gt.GaussianInputFile(file_directory=os.path.join(current_dir, "opt"))

for mol in mols:
    # AllChem.UFFOptimizeMolecule(mol) # not suitable for organometalic complex
    mol = xtb_optimize_molecule(mol)
    g_opt.molecule_to_gjf(mol, cmd=cmd)