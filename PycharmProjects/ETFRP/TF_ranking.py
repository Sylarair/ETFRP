import os
import numpy as np
import logging
import math
import subprocess
import argparse
import sys
from bin import TF_finding_with_peaks, TF_finding_with_reads

logging.basicConfig(level=logging.INFO)


def getargs():
	parser = argparse.ArgumentParser(description='Welcome to use ETFRP to find your TF!')

	current_dir = os.path.realpath(__file__)
	parser.add_argument('-m','--mode', default='peak', action='store', choices=['peak','read'], required=True, dest='')
	parser.add_argument('-i', '--input', required=True, dest='Input file')
	parser.add_argument('-S', '--Signal', required=True, dest='Signal')
	# parser.add_argument('-t', '--treat', dest='treat bam file, can be seperated by comma')
	parser.add_argument('-c', '--control', dest='control bam file, can be seperated by comma')
	parser.add_argument('-o', '--outdir', default=current_dir, dest='Output directory')
	parser.add_argument('-n', '--name', default='Result', dest='Output file name')
	parser.add_argument('-s', '--species', default='hs', choices=['hs', 'mm'], dest='Species')
	parser.add_argument('-r', '--reference', default='hg38', choices=['hg38', 'hg19', 'mm10', 'mm9'], dest='Species')
	# parser.add_argument('-C', '--has_CNV', default=False, type=int, dest='Whether users provide a CNV file')
	parser.add_argument('-C', '--CNV', dest='CNV file name')
	parser.add_argument('-D', '--DEG', dest='Differential expresion genes file name')
	parser.add_argument('-R', '--ReadCount', dest='Genes read counts file name')
	parser.add_argument('-N', '--Top_peak_num', default=2000, type=int, dest='top peak number')

	args = parser.parse_args()

	return args


def load_file(infile):
	if os.file.exists(infile):
		data = [i.strip().split('\t') for i in open(infile, 'r')]
		return data
	else:
		logging.error('{0} does not exist! Please check your file!'.format(infile))

	return True


def bash(*args):
	subprocess.call(*args, shell=True)
	return True


def main():
	args = getargs()
	logging.INFO('Job is starting!')

	top_peak_num = args.Top_peak_num
	if args.input and args.Signal:
		if (args.species == 'hs' and (args.reference == 'hg38' or args.reference == 'hg19')) or (args.species == 'mm' and (args.reference == 'mm10' or args.reference == 'mm9')):
			if args.mode == 'peak':
				# TF_finding_with_peaks(args)
				TF_finding_with_peaks.ranking(args.input, args.outdir, args.name, args.species, args.reference, args.Signal, top_peak_num, args.CNV, args.DEG, args.ReadCount)
			else:
				if args.control:
					logging.INFO('{0} will be considered as control!'.format(args.control))
					TF_finding_with_reads.call_peaks_with_control(args.input, args.control, args.species, args.outdir, args.name)
				else:
					logging.INFO('No control file!'.format(args.control))
					TF_finding_with_reads.call_peaks_without_control(args.input, args.species, args.outdir, args.name)

				peaks = '{0}/{1}_peaks.narrowPeak'.format(args.outdir, args.outprefix)
				TF_finding_with_peaks.ranking(peaks, args.outdir, args.name, args.species, args.reference, args.Signal, top_peak_num, args.CNV, args.DEG, args.ReadCount)

			logging.INFO('Job has been done! Thanks for your using! Bye!')
		else:
			logging.INFO('Your species and reference are not matched! Please check it!')
	else:
		logging.INFO('{0} or {1} is missing! Please check its path!'.format(args.input, args.Signal))

	sys.exit(0)

