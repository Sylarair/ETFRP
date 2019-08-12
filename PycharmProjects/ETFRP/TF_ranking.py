import os
import logging
import subprocess
import argparse
import sys
from bin import TF_finding_with_peaks, TF_finding_with_reads

logging.basicConfig(level=logging.INFO)


def getargs():
	parser = argparse.ArgumentParser(description='Welcome to use ETFRP to find your TF!')

	current_dir = os.path.realpath(__file__)
	parser.add_argument('-m','--mode', default='peak', choices=['peak','read'], required=True, dest='')
	parser.add_argument('-i', '--input', required=True, dest='Input file')
	parser.add_argument('-S', '--Signal', required=True, dest='Signal')
	# parser.add_argument('-t', '--treat', dest='treat bam file, can be seperated by comma')
	parser.add_argument('-c', '--control', action='store_false', dest='control bam file, can be seperated by comma')
	parser.add_argument('-o', '--outdir', default=current_dir, dest='Output directory')
	parser.add_argument('-n', '--name', default='Result', dest='Output file name')
	parser.add_argument('-s', '--species', default='hs', choices=['hs', 'mm'], dest='Species')
	parser.add_argument('-r', '--reference', default='hg38', choices=['hg38', 'hg19', 'mm10', 'mm9'], dest='Species')
	# parser.add_argument('-C', '--has_CNV', default=False, type=int, dest='Whether users provide a CNV file')
	parser.add_argument('-C', '--CNV', action='store_false', dest='CNV file name')
	parser.add_argument('-D', '--DEG', action='store_false', dest='Differential expresion genes files')
	parser.add_argument('-RT', '--TreatReadCount', action='store_false', dest='Genes read counts treat files, will only be used when the DEG file is not specified')
	parser.add_argument('-RC', '--ControlReadCount', action='store_false', dest='Genes read counts control files, will only be used when it is specified, otherwise the control files in database will be used')
	parser.add_argument('-TN', '--TreatName', action='store_false', dest='Treat name is used to determine whether the treat sample is in the control database, will only be used when it is specified, for example: mESC, hindbrain_14.5_day. liver_0_day')
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
	if len(sys.argv) < 4:
		args.print_help(sys.stderr)
		logging.error('Job existing abnormally!')
		sys.exit(1)
	else:
		if args.input and args.Signal:
			if (args.species == 'hs' and (args.reference == 'hg38' or args.reference == 'hg19')) or (args.species == 'mm' and (args.reference == 'mm10' or args.reference == 'mm9')):
				if args.mode == 'peak':
					# TF_finding_with_peaks(args)
					if args.CNV or args.DEG or args.TreatReadCount:
						TF_finding_with_peaks.ranking(args.input, args.outdir, args.name, args.species, args.reference, args.Signal, top_peak_num, args.CNV, args.DEG, args.TreatReadCount, args.ControlReadCount, args.TreatName)
					else:
						logging.error('No CNV or differential expression data or read count data! Please specify at least one of them!')
				else:
					if args.control:
						logging.INFO('{0} will be considered as control!'.format(args.control))
						TF_finding_with_reads.call_peaks_with_control(args.input, args.control, args.species, args.outdir, args.name)
					else:
						logging.INFO('No control file!'.format(args.control))
						TF_finding_with_reads.call_peaks_without_control(args.input, args.species, args.outdir, args.name)

					peaks = '{0}/{1}_peaks.narrowPeak'.format(args.outdir, args.outprefix)
					if args.CNV or args.DEG or args.TreatReadCount:
						TF_finding_with_peaks.ranking(peaks, args.outdir, args.name, args.species, args.reference, args.Signal, top_peak_num, args.CNV, args.DEG, args.TreatReadCount, args.ControlReadCount, args.TreatName)
					else:
						logging.error(
							'No CNV or differential expression data or read count data! Please specify at least one of them!')

				logging.INFO('Job has been done! Thanks for your using! Bye!')
			else:
				logging.INFO('Your species and reference are not matched! Please check it!')
		else:
			logging.INFO('{0} or {1} is missing! Please check its path!'.format(args.input, args.Signal))

		sys.exit(0)

