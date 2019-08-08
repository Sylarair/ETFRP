import os
import numpy as np
import logging
import math
import subprocess
import sys

logging.basicConfig(level=logging.INFO)


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

def call_peaks_with_control(treat_file, control_file, species, outdir, outprefix):
	treat_postfix = treat_file.strip().split('.')[-1]
	control_postfix = control_file.strip().split('.')[-1]

	if treat_postfix == control_postfix and treat_postfix == 'bam':
		macs_mode = 'BAM'
	elif treat_postfix == control_postfix and treat_postfix == 'sam':
		macs_mode = 'SAM'
	elif treat_postfix == control_postfix and treat_postfix == 'BAMPE':
		macs_mode = 'BAMPE'
	elif treat_postfix == control_postfix and treat_postfix == 'BEDPE':
		macs_mode = 'BEDPE'
	else:
		logging.error('Format cannot be processed! Please change another format and try again!')

	bash('macs2 callpeak -t {0} -c {1} -f {2} -g {3} -n {4} --outdir {5} -B --SPMR --nomodel --extsize 147 --shift 73'.format(treat_file, control_file, macs_mode, species, outprefix, outdir))
	# bash('cd {0}; mv {1}_peaks.narrowPeak {1}.bed'.format(outdir, outprefix))



def call_peaks_without_control(treat_file, species, outdir, outprefix):
	treat_postfix = treat_file.strip().split('.')[-1]

	if treat_postfix == 'bam':
		macs_mode = 'BAM'
	elif treat_postfix == 'sam':
		macs_mode = 'SAM'
	elif treat_postfix == 'BAMPE':
		macs_mode = 'BAMPE'
	elif treat_postfix == 'BEDPE':
		macs_mode = 'BEDPE'
	else:
		logging.error('Format cannot be processed! Please change another format and try again!')

	bash('macs2 callpeak -t {0} -f {1} -g {2} -n {3} --outdir {4} -B --SPMR --nomodel --extsize 147 --shift 73'.format(treat_file, macs_mode, species, outprefix, outdir))
	# bash('cd {0}; mv {1}_peaks.narrowPeak {1}.bed'.format(outdir, outprefix))