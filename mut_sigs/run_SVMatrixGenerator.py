#!/usr/bin/env python3

import argparse
from SigProfilerMatrixGenerator.scripts import SVMatrixGenerator as sv

parser = argparse.ArgumentParser(add_help=False)

parser.add_argument('--input_dir', default=None, help='Path to the input folder containing SV files.')
parser.add_argument('--project', default="", help='Name of the project.')
parser.add_argument('--output_dir', default=None, help='Path to the output folder.')

args = parser.parse_args()

sv.generateSVMatrix(input_dir=args.input_dir, project=args.project, output_dir=args.output_dir)