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


def download_reference(reference, outpath):
	if not os.path.exists(reference):
		reference_name = reference.strip().split('/')[-1].split('.')[0]
		bash('mkdir {0}; wget -P {0} -c https://hgdownload.soe.ucsc.edu/goldenPath/{1}/bigZips/{1}.fa.gz; cd {0}; gunzip {1}.fa.gz'.format(outpath, reference_name))


def rank_and_select_peaks(inpeaks, signal, top_peak_num, outdir, outprefix):
	# peaks_signal_rank = [[i[0], i[1], i[2], str(i[0]) + str(i[1]) + str(i[2])] for i in peaks]
	peaks_signal_rank = '{0}/{1}.bed'.format(outdir, outprefix)
	bash("""awk '{{print $1"\t"$2"\t"$3"\t"NR}}' {0} > {1}""".format(inpeaks, peaks_signal_rank))
	bash('sort -k4,4 {0} > {0}.temp'.format(peaks_signal_rank))
	bash('bigWigAverageOverBed {0} {1}.temp {1}.temp2'.format(signal, peaks_signal_rank))
	bash('sort -1,1 {0}.temp2 | cut -f4 > {0}.temp3; paste {1} {0}.temp3 | sort -k5gr,5gr > {0}.final_sorted.bed'.format(peaks_signal_rank))
	bash('head -n {0} {1}.final_sorted.bed > {1}.final_sorted.{0}.bed'.format(top_peak_num,
			peaks_signal_rank))


def AME_enrich_TFs(inpeaks, outdir, outprefix, top_peak_num, species, reference):

	current_dir = os.path.realpath(__file__)
	reference_file = '{0}/library/{1}/{2}.fa'.format(current_dir, species, reference)
	PWM_file = '{0}/library/{1}/{1}.total_motif_PWM.meme'.format(current_dir, species)
	outpath = '{0}/{1}_AME_results'.format(outdir, outprefix)
	bash('mkdir -p {0}'.format(outpath))
	# top_peaks = '{0}/{1}.final_sorted.top{2}.bed'.format(outdir, outprefix, top_peak_num)
	top_peaks_fa = '{0}/{1}.final_sorted.top{2}.fa'.format(outdir, outprefix, top_peak_num)
	background_peaks_fa = '{0}/{1}.background.fa'.format(outdir, outprefix)
	background_info = '{0}/{1}.background.info'.format(outdir, outprefix)
	control_fa = '{0}/{1}.control.fa'.format(outdir, outprefix)
	bash('head -n {0} | bedtools getfasta -fi {1} -bed stdin -fo {2}'.format(inpeaks, reference_file, top_peaks_fa))
	bash('tail -n +{0} {1} | tail -n 10000 | bedtools getfasta -fi {2} -bed stdin -fo {3}'.format(int(top_peak_num) + 1, inpeaks, reference_file, background_peaks_fa))
	bash('fasta-shuffle-letters {0} -dna -seed 1 {1}'.format(background_peaks_fa, control_fa))
	bash('fasta-get-markov {0} -dna -m 0 {1}'.format(background_peaks_fa, background_info))

	bash('ame --oc ${0} --control {1} --bfile {2} {3} {4}'.format(outpath, control_fa, background_peaks_fa, top_peaks_fa, PWM_file))
	bash('bash choose_best_motif.sh {0} {1}'.format(outpath, outdir)) ***


def AME_pvalue_CNV_rank_product(AME_enriched_TFs3, CNV_file, out_TFs_file):
	CNV = load_file(CNV_file)  # CNV file must be in a format like: gene1 (hugo id) CNV_ratio
	# CNV2 = [[CNV[i][0], CNV[i][1], i + 1] for i in range(len(CNV))]
	CNV_dict = {}
	for i in range(len(CNV)):
		if CNV[i][0] not in CNV_dict.keys():
			CNV_dict[CNV[i][0]] = [CNV[i][1], i + 1]

	merge_list = []
	for i in range(len(AME_enriched_TFs3)):
		_item = AME_enriched_TFs3[i]
		TF = _item[0]
		if TF in CNV_dict.keys():
			merge_list.append([TF, _item[1], _item[2], CNV_dict[TF][0], CNV_dict[TF][1]])
		else:
			merge_list.append([TF, _item[1], _item[2], np.nan, np.nan])

	merge_list2 = []
	max_CNV_rank = np.nanmax([i[4] for i in merge_list])

	for i in range(len(merge_list)):
		CNV_rank = merge_list[i][4]
		if CNV_rank != np.nan:
			merge_list2.append(merge_list[i] + [merge_list[i][2] * merge_list[i][4]])
		else:
			merge_list2.append(merge_list[i] + [merge_list[i][2] * (max_CNV_rank + 1)])

	merge_list3 = sorted(merge_list2, lambda x: x[5])

	with open(out_TFs_file, 'w') as g:
		for _line in merge_list3:
			print >> g, '\t'.join(np.array(_line, dtype=str).tolist())


def run_DESeq2(read_count_file_folder, DEG_file):
	***


def AME_pvalue_TF_gene_FC_rank_product(AME_enriched_TFs3, DEG_file, read_count_file_folder, out_TFs_file):
	if not os.path.exists(DEG_file):
		run_DESeq2(read_count_file_folder, DEG_file)

	DEG = load_file(DEG_file)  # DEG file must be in a format like: gene1 (hugo id) baseMean log2FoldChange lfcSE pvalue padj
	# DEG2 = [[DEG[i][0], DEG[i][1], i + 1] for i in range(len(DEG))]
	DEG_dict = {}
	for i in range(len(DEG)):
		if DEG[i][0] not in DEG_dict.keys():
			DEG_dict[DEG[i][0]] = [DEG[i][2], i + 1]

	merge_list = []
	for i in range(len(AME_enriched_TFs3)):
		_item = AME_enriched_TFs3[i]
		TF = _item[0]
		if TF in DEG_dict.keys():
			merge_list.append([TF, _item[1], _item[2], DEG_dict[TF][0], DEG_dict[TF][1]])
		else:
			merge_list.append([TF, _item[1], _item[2], np.nan, np.nan])

	merge_list2 = []
	max_DEG_rank = np.nanmax([i[4] for i in merge_list])

	for i in range(len(merge_list)):
		DEG_rank = merge_list[i][4]
		if DEG_rank != np.nan:
			merge_list2.append(merge_list[i] + [merge_list[i][2] * merge_list[i][4]])
		else:
			merge_list2.append(merge_list[i] + [merge_list[i][2] * (max_DEG_rank + 1)])

	merge_list3 = sorted(merge_list2, lambda x: x[5])

	with open(out_TFs_file, 'w') as g:
		for _line in merge_list3:
			print >> g, '\t'.join(np.array(_line, dtype=str).tolist())


def ranking(infile, outdir, outprefix, species, reference, signal, top_peak_num, CNV_file, DEG_file, read_count_file_folder):
	peaks = load_file(infile)
	current_dir = os.path.realpath(__file__)
	reference = '{0}/library/{1}/{2}.fa'.format(current_dir, species, reference)
	reference_path = '{0}/library/{1}'.format(current_dir, species)

	download_reference(reference, reference_path)

	# rank by signal and select top 2000 peaks
	rank_and_select_peaks(peaks, signal, top_peak_num, outdir, outprefix)


	ranked_peaks = '{0}/{1}.final_sorted.bed'.format(outdir, outprefix, top_peak_num)
	AME_enrich_TFs(ranked_peaks, outdir, outprefix)

	AME_enriched_TFs_file = '{0}/{1}.total_motif.txt'.format(outdir, outprefix)
	AME_enriched_TFs0 = load_file(AME_enriched_TFs_file)
	AME_enriched_TFs = [[i[3], float(i[5])] for i in AME_enriched_TFs0] # col4,6 are TFs id and AME_pvalue
	AME_enriched_TFs2 = sorted(AME_enriched_TFs, lambda x: x[1])
	AME_enriched_TFs3 = [[AME_enriched_TFs2[i][0],AME_enriched_TFs2[i][1], i+1] for i in range(len(AME_enriched_TFs2))] # add rank

	out_TFs_file = '{0}/{1}.found_TFs.txt'.format(outdir, outprefix)
	if CNV_file:
		AME_pvalue_CNV_rank_product(AME_enriched_TFs3, CNV_file, out_TFs_file)
	else:
		AME_pvalue_TF_gene_FC_rank_product(AME_enriched_TFs3, DEG_file, read_count_file_folder, out_TFs_file)


