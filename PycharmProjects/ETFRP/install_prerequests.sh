local_path=`pwd`
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/ ${local_path}/bin/shell_scripts
export PATH:${local_path}/bin/shell_scripts:$PATH
conda install -f requirements.txt