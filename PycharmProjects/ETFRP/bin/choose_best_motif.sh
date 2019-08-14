#!/bin/bash
inpath=$1
outpath=$2
prefix=$3
species=$4

function choose_most_significant_one_motif_for_TF(){
	inpath=$1
	outpath=$2
	prefix=$3
	species=$4

	outprefix=${outpath}"/"${prefix}
	#
	total_motif_temp=${outprefix}".total-motif.temp.bed"
	total_motif=${outprefix}".total-motif.bed"
	[[ -s ${total_motif} ]] && rm ${total_motif}

	cat ${inpath}"/ame.tsv" | grep -v 'motif_DB' | grep -v '^#' | grep -v '^\s*$' | sort -k7g,7g > ${total_motif_temp}
	for k in `cut -f4 ${total_motif_temp} | sort | uniq`
	do
		awk -v motif=${k} '{if($4==motif) print $0}' ${total_motif_temp} | sort -k7g,7g > ${total_motif}".temp"
		head -n1 ${total_motif}".temp" >> ${total_motif}
	done

	sort -k7g,7g ${total_motif} > ${total_motif}".temp2"
	mv ${total_motif}".temp2" ${total_motif}
	rm -f ${total_motif_temp}

	#
	if [[ ${species} == "hs" ]];then
	  sed 's/(//g' ${total_motif} | sed 's/)//g' | awk '{split($4,a,"_");print toupper(a[1])"\t"$5"\t"$6"\t"$7"\t"$8"\t"$2}' > ${outprefix}".total-motif.main-info.temp.bed"
	else
	  sed 's/(//g' ${total_motif} | sed 's/)//g' | awk '{split($4,a,"_");print toupper(a[1])"\t"$5"\t"$6"\t"$7"\t"$8"\t"$2}' | awk '{printf "%s\t%s\t";toupper(substr($1,1,1)),tolower(substr($1,2));print $2"\t"$3"\t"$4"\t"$5"\t"$6}' > ${outprefix}".total-motif.main-info.temp.bed"
	fi

	summary_temp_file=${outprefix}".total-motif.main-info.temp.bed"
	summary_file=${outprefix}".total-motif.main-info.bed"
	[[ -s ${summary_file} ]] && rm -f ${summary_file}
	for i in `cut -f1 ${summary_temp_file} | sort | uniq`
	do
		awk -v motif=${i} '{if($1==motif) print $0}' ${summary_temp_file} | sort -k3g,3g > ${summary_temp_file}".temp"
		head -n1 ${summary_temp_file}".temp" >> ${summary_file}
	done

	sort -k3g,3g ${summary_file} > ${summary_file}".sorted"
	cat ${summary_file}".sorted" > ${summary_file}
}

choose_most_significant_one_motif_for_TF ${inpath} ${outpath} ${prefix} ${species}