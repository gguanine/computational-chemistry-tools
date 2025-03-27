import cc_tools.gaussian_tools as gt
import os

current_dir = os.path.dirname(os.path.abspath(__file__))

input_dir = os.path.join(current_dir, "opt") # input directory for reading Gaussian output files.
output_dir = os.path.join(current_dir, "energy") # output directory for Gaussian calculation.
cmd = "# m062x/6-311g** SCRF(SMD, Solvent=Dichloromethane)" # command for energy computation

if __name__ == "__main__":
    gt.dir_logs_to_gjf(input_dir, output_dir, cmd, name_suffix='opt')