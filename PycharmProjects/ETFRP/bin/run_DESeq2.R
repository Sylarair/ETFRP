if (! require('DESeq2')){
    if (! require('BiocManager')){
        install.packages('BiocManager')
        suppressMessages(library('BiocManager'))
        BiocManager::install('DESeq2')
    }
}
if (! require('dplyr')){
    install.packages('dplyr')
}

if (! require('tximport')){
    BiocManager::install('tximport')
}

suppressMessages(library('DESeq2'))
suppressMessages(library('dplyr'))
suppressMessages(library('tximport'))

args = commandArgs(T)
current_dir = args[1]
read_count_folder = args[2]
results = args[3]
species = args[4]
control_read_count_folder = args[5]
treat_name = args[6]

check_is_readcount_file = function(infolder){
    files = list.files(infolder)
    random_choosen_file = paste0(infolder,'/',files[sample(c(1:length(files)),1)])
    print(random_choosen_file)
    random_choosen_file_concept = read.table(random_choosen_file, header=T, sep='\t', stringsAsFactors=F)
    random_choosen_file_concept[,1] = as.character(as.vector(sapply(as.vector(random_choosen_file_concept[,1],function(x) strsplit(x,'\\.')[[1]][1])))) # remove Ensembl id version number, eg. ENSG00000000003.10 --> ENSG00000000003
    random_choosen_file_concept = random_choosen_file_concept[order(-random_choosen_file_concept[,2]),] # sort data frame by the read count number from the largest to smallest

    if (random_choosen_file_concept[1,2] %% 1 == 0 & random_choosen_file_concept[sample(1,nrow(random_choosen_file_concept),1),2] %% 1 == 0){
        return('TRUE')
    } else {
        return('FALSE')
    }
}

check_is_readcount_file2 = function(infolder) {
  files = infolder #list.files(infolder)
  random_choosen_file = files[sample(c(1:length(files)), 1)] #paste0(infolder, '/', files[sample(1, length(files), 1)])
  print(random_choosen_file)
  random_choosen_file_concept = read.table(random_choosen_file, header = T, sep = '\t', stringsAsFactors = F)
  random_choosen_file_concept[, 1] = as.character(as.vector(sapply(as.vector(random_choosen_file_concept[, 1]), function(x) strsplit(x, '\\.')[[1]][1]))) # remove Ensembl id version number, eg. ENSG00000000003.10 --> ENSG00000000003
#   print(str(random_choosen_file_concept))
  random_choosen_file_concept = random_choosen_file_concept[order(-random_choosen_file_concept[, 6]),] # sort data frame by the read count number from the largest to smallest

  if (random_choosen_file_concept[1, 6] %% 1 == 0 & random_choosen_file_concept[sample(c(1:nrow(random_choosen_file_concept)), 1), 6] %% 1 == 0) {
    return('TRUE')
  } else {
    return('FALSE')
  }
}

compute_DEG = function(current_dir, infolder, species, output){
    # files = list.files(infolder)
    treat_files = as.character(as.vector(strsplit(infolder,',')[[1]])) #paste0(infolder,'/',files)
    print(infolder)
    print(treat_files)
    if (control_read_count_folder=='False'){
        control_files = paste0(current_dir,'/library/Tissue_read_counts/',species,'/',list.files(paste0(current_dir,'/library/Tissue_read_counts/', species)))
    } else{
        control_files = paste0(control_read_count_folder, '/', list.files(control_read_count_folder)) # actually it is control files
    }
    # print(control_files)
    if (treat_name!='False' & control_read_count_folder=='False'){
    #   print(grep(treat_name, control_files))
       control_files = control_files[grep(treat_name,control_files, invert = TRUE)]
    }
    # print(control_files)
    # print(treat_name)
    # print(length(treat_files))
    # print(length(control_files))
    if (length(treat_files) != 0 & length(control_files) != 0){
        # load tissue / cell lines rsem results and convert them into read counts
        if (check_is_readcount_file2(control_files) == 'True'){
            sample_files = control_files #unlist(c(treat_files, control_files))
            sample_condition = rep('control',length(control_files)) #c(rep('treat',length(treat_files)), rep('control',length(control_files)))

            sampleTable = data.frame(sampleName = as.character(as.vector(sapply(sample_files,function(x) basename(x)))), fileName = as.character(as.vector(sapply(sample_files,function(x) basename(x)))), condition = sample_condition)
            ddsHTSeq = DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = dirname(control_files[1]), design= ~ condition)
            control_counts = assay(ddsHTSeq)
        } else{
            txrsem = tximport(control_files,type='rsem',txIn=FALSE,txOut=F)
            row_num <- which(apply(txrsem$length,1,min)==0)
            txrsem$length <- txrsem$length[-1*row_num,]
            txrsem$abundance <- txrsem$abundance[-1*row_num,]
            txrsem$counts <- txrsem$counts[-1*row_num,]
            # print(str(txrsem))
            # print(str(txrsem$counts))
            # print(str(rownames(txrsem$counts)))
            # print(str(names(txrsem$counts)))
            # print(str(rownames(txrsem$counts)))
            control_counts = as.data.frame(txrsem$counts)
        }
        # load treat files
        if (check_is_readcount_file2(infolder) == 'TRUE'){
            sample_files = treat_files #unlist(c(treat_files, control_files))
            sample_condition = rep('treat',length(treat_files)) #c(rep('treat',length(treat_files)), rep('control',length(control_files)))

            sampleTable = data.frame(sampleName = as.character(as.vector(sapply(sample_files,function(x) basename(x)))), fileName = as.character(as.vector(sapply(sample_files,function(x) basename(x)))), condition = sample_condition)
            ddsHTSeq = DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = dirname(sample_files[1]), design= ~ condition)
            treat_counts = as.data.frame(assay(ddsHTSeq))
            # print(str(treat_counts))
        } else {
            treat_txrsem = tximport(treat_files,type='rsem',txIn=FALSE,txOut=F)
            treat_row_num <- which(apply(treat_txrsem$length,1,min)==0)
            treat_txrsem$length <- treat_txrsem$length[-1*treat_row_num,]
            treat_txrsem$abundance <- treat_txrsem$abundance[-1*treat_row_num,]
            treat_txrsem$counts <- treat_txrsem$counts[-1*treat_row_num,]
            treat_counts = as.data.frame(treat_txrsem$counts)
            # print(str(rownames(treat_txrsem$counts)))
            # print(str(names(treat_txrsem$counts)))
            if (ncol(treat_counts) > 1){
                rownames(treat_counts) = rownames(treat_txrsem$counts)
            } else{
                rownames(treat_counts) = names(treat_txrsem$counts)
            }
            # print(str(treat_counts))
        }

        control_counts$gene_id = rownames(control_counts)
        treat_counts$gene_id = rownames(treat_counts)
        # print(str(treat_counts))
        # print(str(control_counts))
        total_counts = full_join(treat_counts, control_counts, by='gene_id')
        total_counts = unique(total_counts[!is.na(apply(total_counts,1,max)),])
        rownames(total_counts) = total_counts$gene_id
        total_counts = total_counts[, colnames(total_counts)!='gene_id']
        colnames(total_counts) = c(paste0('treat', c(1:(ncol(treat_counts) - 1))), paste0('control', c(1:(ncol(control_counts) - 1))))
        total_counts = round(total_counts)
        print(str(total_counts))
        colData = data.frame(row.names = colnames(total_counts), condition = unlist(c(rep('treat',ncol(treat_counts)-1),rep('control',ncol(control_counts)-1))))
        dds = DESeqDataSetFromMatrix(countData = total_counts, colData = colData, design = ~ condition)
        dds$condition = factor(dds$condition, levels = c('control','treat')) # set levels to make treat compared to control
        dds = DESeq(dds)
        res = results(dds)
        # print(str(res))
        res2 = as.data.frame(res@listData)
        rownames(res2) = res@rownames
        res = res2
        res = res[!is.na(res$log2FoldChange),]
        # print(str(res))
        #res_mat = data.frame(gene_id=rownames(res),fc=res[,2],adjust_pvalue=res[,6],pvalue=res[,5])
        res_mat = res
        res_mat[,5] = res_mat[,5]+10^(-300) # avoid some zero
        res_mat[,6] = res_mat[,6]+10^(-300) # avoid some zero
        res_mat = res_mat[order(-res_mat$log2FoldChange, res_mat$pvalue),]
        res_mat$ensembl_id = as.character(as.vector(sapply(as.vector(rownames(res_mat)), function(x) strsplit(x,'\\.')[[1]][1])))
        # print(str(res_mat))
        id_transversion = read.table(paste0(current_dir,'/library/',species,'.gene_info.txt'), sep='\t', header=F, stringsAsFactors=F)
        colnames(id_transversion) = c('chromosome', 'source', 'type', 'start', 'end', 'strand', 'ensembl_id_with_version', 'hugo_id', 'ensembl_id')
        res_mat = left_join(res_mat, id_transversion, by='ensembl_id')
        # print(str(res_mat))
        res_mat = res_mat[,c(15,c(1:7))]
        res_mat = res_mat[,-5]
        write.table(res_mat, file=output, sep='\t', row.names=F, col.names=T, quote=F)
    } else{
        print('No treat files or control files! Please check it!')
    }
}

compute_DEG(current_dir, read_count_folder, species, results)