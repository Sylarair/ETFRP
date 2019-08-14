import os
import logging
import subprocess
import argparse
import sys
from bin import TF_finding_with_peaks, TF_finding_with_reads

logging.basicConfig(format='%(asctime)s.%(msecs)03d %(levelname)s {%(module)s} [%(funcName)s] %(message)s',
                    datefmt='%Y-%m-%d,%H:%M:%S', level=logging.INFO)
logger = logging.getLogger()
current_dir = os.path.realpath(__file__)


# print(current_dir)

# def getargs():
# 	parser = argparse.ArgumentParser(description='Welcome to use ETFRP to find your TF!')
#
# 	parser.add_argument('-m','--mode', default='peak', choices=['peak','read'], required=True, dest='')
# 	parser.add_argument('-i', '--input', required=True, dest='Input file')
# 	parser.add_argument('-S', '--Signal', required=True, dest='Signal')
# 	# parser.add_argument('-t', '--treat', dest='treat bam file, can be seperated by comma')
# 	parser.add_argument('-c', '--control', default=False,
# 	                    dest='control bam file, can be seperated by comma')
# 	parser.add_argument('-o', '--outdir', default=current_dir, dest='Output directory')
# 	parser.add_argument('-n', '--name', default='Result', dest='Output file name')
# 	parser.add_argument('-s', '--species', default='hs', choices=['hs', 'mm'], dest='Species')
# 	parser.add_argument('-r', '--reference', default='hg38', choices=['hg38', 'hg19', 'mm10', 'mm9'], dest='Species')
# 	parser.add_argument('-N', '--Top_peak_num', nargs=1, type=int, default=2000, dest='top peak number')
# 	# parser.add_argument('-C', '--has_CNV', default=False, type=int, dest='Whether users provide a CNV file')
# 	parser.add_argument('-C', '--CNV', default='False', dest='CNV file name')
# 	parser.add_argument('-D', '--DEG', default='False', dest='Differential expresion genes files')
# 	parser.add_argument('-RT', '--TreatReadCount', default='False', dest='Genes read counts treat files, will only be used when the DEG file is not specified')
# 	parser.add_argument('-RC', '--ControlReadCount', default='False', dest='Genes read counts control files, will only be used when it is specified, otherwise the control files in database will be used')
# 	parser.add_argument('-TN', '--TreatName', default='False', dest='Treat name is used to determine whether the treat sample is in the control database, will only be used when it is specified, for example: mESC, hindbrain_14.5_day. liver_0_day')
#
# 	args = parser.parse_args()
#
# 	return args


def getargs2():
	parser = argparse.ArgumentParser(description='Welcome to use ETFRP to find your TF!')
	parser.add_argument('-m', dest='mode', default='peak', choices=['peak', 'read'], required=True, help='')
	parser.add_argument('-i', dest='input', required=True, help='Input file')
	parser.add_argument('-S', dest='Signal', required=True, help='Signal')
	# parser.add_argument('-t', '--treat', dest='treat bam file, can be seperated by comma')
	parser.add_argument('-c', dest='control', nargs='+', default=False, help='control bam file, can be seperated by comma')
	parser.add_argument('-o', dest='outdir', nargs="?", default=current_dir, help='Output directory')
	parser.add_argument('-n', dest='name', nargs="?", default='Result', help='Output file name')
	parser.add_argument('-s', dest='species', nargs="?", default='hs', choices=['hs', 'mm'], help='Species')
	parser.add_argument('-r', dest='reference', nargs="?", default='hg38', choices=['hg38', 'hg19', 'mm10', 'mm9'], help='reference')
	parser.add_argument('-N', dest='Top_peak_num', nargs="?", type=int, default=2000, help='top peak number')
	parser.add_argument('-C', dest='CNV', nargs="?", default='False', help='CNV file name')
	parser.add_argument('-D', dest='DEG', nargs="?", default='False', help='Differential expresion genes files')
	parser.add_argument('-RT', dest='TreatReadCount', nargs="?", default='False', help='Genes read counts treat files, will only be used when the DEG file is not specified')
	parser.add_argument('-RC', dest='ControlReadCount', nargs="?", default='False', help='Genes read counts control files, will only be used when it is specified, otherwise the control files in database will be used')
	parser.add_argument('-TN', dest='TreatName', nargs="?", default='False', help='Treat name is used to determine whether the treat sample is in the control database, will only be used when it is specified, for example: mESC, hindbrain_14.5_day. liver_0_day')

	args = parser.parse_args()

	return args


def load_file(infile):
	if os.path.exists(infile):
		data = [i.strip().split('\t') for i in open(infile, 'r')]
		return data
	else:
		logger.error('{0} does not exist! Please check your file!'.format(infile))

	return True


def bash(args):
	# subprocess.call(*args, shell=True)
	# return True
	return os.system(args)


def main():
	args = getargs2()
	logger.info('Job is starting!')
	print(sys.argv)
	top_peak_num = args.Top_peak_num
	if len(sys.argv) < 4:
		# args.print_help(sys.stderr)
		args.print_help()
		logging.error('Job existing abnormally!')
		sys.exit(1)
	else:
		if args.input != 'False' and args.Signal != 'False':
			if (args.species == 'hs' and (args.reference == 'hg38' or args.reference == 'hg19')) or (
					args.species == 'mm' and (args.reference == 'mm10' or args.reference == 'mm9')):
				if args.mode == 'peak':
					# TF_finding_with_peaks(args)
					if args.CNV != 'False' or args.DEG != 'False' or args.TreatReadCount != 'False':
						TF_finding_with_peaks.ranking(args.input, args.outdir, args.name, args.species, args.reference, args.Signal, top_peak_num, args.CNV, args.DEG, args.TreatReadCount, args.ControlReadCount, args.TreatName)
					else:
						logger.error(
							'No CNV or differential expression data or read count data! Please specify at least one of them!')
				else:
					if args.control != 'False':
						logger.info('{0} will be considered as control!'.format(args.control))
						TF_finding_with_reads.call_peaks_with_control(args.input, args.control, args.species, args.outdir, args.name)
					else:
						logger.info('No control file!'.format(args.control))
						TF_finding_with_reads.call_peaks_without_control(args.input, args.species, args.outdir, args.name)

					peaks = '{0}/{1}_peaks.narrowPeak'.format(args.outdir, args.name)
					if args.CNV != 'False' or args.DEG != 'False' or args.TreatReadCount != 'False':
						TF_finding_with_peaks.ranking(peaks, args.outdir, args.name, args.species,args.reference, args.Signal, top_peak_num, args.CNV, args.DEG, args.TreatReadCount, args.ControlReadCount, args.TreatName)
					else:
						logger.error('No CNV or differential expression data or read count data! Please specify at least one of them!')

				logger.info('Job has been done! Thanks for your using! Bye ;)!')
			else:
				logger.info('Your species and reference are not matched! Please check it!')
		else:
			logger.info('{0} or {1} is missing! Please check its path!'.format(
				args.input, args.Signal))

		sys.exit(0)


if __name__ == '__main__':
	main()
