import cc_tools.gaussian_tools as gt
import argparse
import os

if __name__ == "__main__":
    current_dir = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--input", required=True, help="input directory")
    parser.add_argument("-o", "--output", required=True, help="output directory")
    parser.add_argument("-c", "--cmd", required=True, help="computation command")
    args = parser.parse_args()
    
    input_dir = os.path.join(current_dir, args.input)
    output_dir = os.path.join(current_dir, args.output)
    
    gt.dir_logs_to_gjf(input_dir, output_dir, args.cmd, name_suffix=args.input)