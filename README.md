# del_score_pipe
This code is used to score the deletion, which is based on the secondary structure of DNA near the genome of the deletion, and GC content.
###file: chr,start,end,strand
Rscript get_bed.R info.csv
###discard "
sed -i 's/"//g' bed_file.txt
###convert to bed
awk 'BEGIN{ FS=" ";OFS="\t" }{ print $1,$2,$3,$4,$5,$6 }' bed_file.txt > bed_file.bed
####Note that some positions are represented by scientific counting

###getfasta
###the reference fasta file can be downloaded at http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz
###bedtools v2.27.1
bedtools getfasta -fi GRCh37.primary_assembly.genome.fa -bed bed_file.bed -s -name -fo del_flank.fa
###non-B gfa can be downloaded at https://nonb-abcc.ncifcrf.gov/apps/site/default
non-B_gfa/gfa -seq del_flank.fa -out del
Rscript score.R info.csv
