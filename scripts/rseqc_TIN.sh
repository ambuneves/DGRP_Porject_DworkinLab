### Script to calculate the transcript integrity number (TIN) using RSeQC --> tin.py

#The files that tin.py uses are a reference file (bed12) and an alignment file (.bam) that has been sorted and indexed (using SAMtools)

######### I've already generated the index [/2/scratch/amandaN/index_files/dmel-all-r6.38_good.bed], so you can skip this but here's how to do it for another version of the transcriptome

#You will need a reference which is in 12 column BED format (.bed_)
#You can download one from FlyBase, it will be a (.gtf) file which you will need to convert

#First download UCSC package
wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
chmod 755 gtfToGenePred genePredToBed

#Convert gtf to GenePred and then convert the GenePred to bed12

gtfToGenePred reference.gtf reference.genePhred #convert from gtf to genePhred
genePredToBed reference.genePhred reference.bed #convert from genePhred to bed12
rm reference.genePhred #remove the intermediate genePhred file

#The reference I am using is: /2/scratch/amandaN/index_files/dmel-all-r6.38_good.bed

#### When running STAR to align, you can set it to ouput a sorted .bam file, it will have the suffix [.sortedByCoord.out.bam]
#This file is sorted, you need to use SAMtools to index it.
#In a screen, index all of the files at once:

ls *sortedByCoord.out.bam | xargs -n1 -P5 samtools index

############# Once you've got the sorted and indexed .bam files, you can run tin.py!
# For a single file:

python3 /usr/local/bin/tin.py -i [file].sortedByCoord.out.bam -r /2/scratch/amandaN/index_files/dmel-all-r6.38_good.bed


# For all of the files in a directory, RSeQC takes as an input a directory with BAM files, so give it the name of your directory:

python3 /usr/local/bin/tin.py -i /path/to/directory -r /2/scratch/amandaN/index_files/dmel-all-r6.38_good.bed
