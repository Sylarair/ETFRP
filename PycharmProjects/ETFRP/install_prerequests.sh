rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/ ./bin/shell_scripts
export PATH:`pwd`/bin/shell_scripts:$PATH
conda install -f requirements.txt