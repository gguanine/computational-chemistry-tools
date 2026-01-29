import subprocess
import shutil
import tempfile
import os
from rdkit import Chem
from . import utils

def xtb_optimize_molecule(mol, cmd:list=[]):
    '''
    use xtb to optimize a molecule. 
    Charge and multiplet are extracted from mol obj.
    Extra command is listed in cmd as a list
    '''
    
    # check if xtb is available
    path_xtb =  shutil.which("xtb")
    if path_xtb:
        print("found xtb binary at location", path_xtb)
    else:
        raise FileNotFoundError("xtb binary not found. Plesae add to PATH or configurate at conig.py")

    # Create a temporary directory for xtb runing
    with tempfile.TemporaryDirectory() as temp_dir:
        input_file = os.path.join(temp_dir, "input.xyz")
        output_file = os.path.join(temp_dir, "xtbopt.xyz")  # If the binary generates output

        # Write some data to an input file
        with open(input_file, "w") as f:
            f.write(Chem.MolToXYZBlock(mol))

        
        command = ["xtb", input_file, "--opt"] + cmd
        name = mol.GetProp("_Name")
        c = mol.GetProp("Charge") # charge
        u = int(mol.GetProp("Multiplicity")) - 1 # unpaired electrons
        
        if c != 0:
            command += ["-c", str(c)]
            
        if u != 0:
            command += ["-u", str(u)]
        
        # Call the external binary and pass file paths
        try:
            subprocess.run(
                command,
                capture_output=True,
                text=True,
                check=True,
                encoding="utf-8",
                cwd=temp_dir
            )

            # Read and process output file if needed
            if os.path.exists(output_file):
                mol2 = Chem.MolFromXYZFile(output_file)
                utils.set_cc_props(mol2, name=name, charge=c, multiplicity=u+1)
                print(f"Sucessfully optimize {name} using xtb")
                return mol2

        except subprocess.CalledProcessError as e:
            print("Error:", e)
            
    print(f"Unable to optimize {name} with xbt")
            
