#!/usr/bin/python

''' 
Script to create scratch files and submit jobs for bowtie2 mapping of ChIP-seq or DamID-seq data.
Description- Creates and submit pbs files for all the .fq files in current directory
converts sam output to bam format, sorts bam files
and removes duplicates using picard. Files with extension .sorted.bam files are original alignment files and
.nodup.bam are after removing duplicates. Log file for this step is log_bowtie.txt.

Copy do_bowtiemapping.py and designfile.txt to current folder
USAGE: python do_bowtiemapping.py -d designfile.txt -p /path/chipseqscripts

chipseqscripts is the folder containing the bowtie indexes and Araport 11 gff file

1) Creates a log file for all the commands that ran for each sample

'''

import sys,os
import shutil
import glob
from subprocess import call
import argparse as agp

def domapping(designfile,scriptpath):

	logfile= open('log_bowtie.txt','wb')

	if scriptpath[-1]=='/':
        	scriptpath= scriptpath[:-1]

        datadir= os.getcwd()
        print 'Current working directory= ',datadir
        for myfile in glob.glob('*.fq'):
                print 'myfile= ',myfile
                sbatch= open(datadir+'/sbatch_bowtiemapping.sh','w')
                sbatch.write('#!/bin/sh\n')
                sbatch.write('#\n')
                sbatch.write('#SBATCH --verbose\n')
                sbatch.write('#SBATCH --job-name=run_bowtie\n')
                sbatch.write('#SBATCH --output=bowtie_%j.out\n')
                sbatch.write('#SBATCH --error=bowtie_%j.err\n')
                sbatch.write('#SBATCH --time=01:00:00\n')
                sbatch.write('#SBATCH --nodes=1\n')
                sbatch.write('#SBATCH --cpus-per-task=4\n')
                sbatch.write('#SBATCH --mem=20GB\n')
                sbatch.write('\n')
                sbatch.write('cd '+datadir+'\n')
                sbatch.write('module load bowtie2/intel/2.2.9\n')
                sbatch.write('module load samtools/intel/1.3.1\n')
		sbatch.write('module load picard-tools/1.88\n')
                mycmd= ''.join(['bowtie2 -p 4 -t -x ',scriptpath,'/bowtie_index_Ath/TAIR10_Ath_Allchr -U ',myfile.strip(),' -S ',myfile.strip().split('.')[0]+'.sam\n'])
                sbatch.write(mycmd)
                sbatch.write('# Convert sam to bam files\n')
		sbatch.write(''.join(['samtools view -bS ',(myfile.strip().split('.')[0]+'.sam'),' > ',(myfile.strip().split('.')[0]+'.bam'),'\n']))
		sbatch.write(''.join(['samtools sort ',(myfile.strip().split('.')[0]+'.bam'),' -o ',(myfile.strip().split('.')[0]+'.sorted.bam'),'\n']))
		sbatch.write(''.join(['samtools index ',(myfile.strip().split('.')[0]+'.sorted.bam'),'\n']))  
		sbatch.write('# Mark duplicates using picard tools\n')
		picardcmd= ''.join(['java -jar /share/apps/picard/2.8.2/picard-2.8.2.jar MarkDuplicates I=',(myfile.strip().split('.')[0]+'.sorted.bam'),' O=',(myfile.strip().split('.')[0]+'.markdup.bam'),' M=',(myfile.strip().split('.')[0]+'.marked.txt'),'\n'])
		sbatch.write(picardcmd)
		sbatch.write('#Remove duplicates\n')
		sbatch.write(''.join(['samtools view -bF 1028 ',(myfile.strip().split('.')[0]+'.markdup.bam'),' > ',(myfile.strip().split('.')[0]+'.nodup.bam'),'\n']))
		sbatch.write(''.join(['samtools index ',(myfile.strip().split('.')[0]+'.nodup.bam'),'\n']))
		sbatch.write('rm '+(myfile.strip().split('.')[0]+'.markdup.bam')+'\n')
		sbatch.write('rm '+(myfile.strip().split('.')[0]+'.marked.txt')+'\n')
                sbatch.write('\n')
		sbatch.close()
                shutil.copyfileobj(open("sbatch_bowtiemapping.sh",'rb'), logfile)
		call(["sbatch", "sbatch_bowtiemapping.sh"])
			
	
if __name__=='__main__':
        parser= agp.ArgumentParser()
        parser.add_argument('-d',help='Provide design file for your data')
        parser.add_argument('-p',help='Provide path for scripts')
        args = parser.parse_args()
        domapping(args.d,args.p)


