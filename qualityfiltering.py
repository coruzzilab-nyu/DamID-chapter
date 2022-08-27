#!/usr/bin/python
''' 
This script is to rename files, to create sbatch files and submit job for quality filtering of ChIP-seq or DamID-seq data.
Description- Creates and submit pbs files for all the .fq files in current directory

Copy qualityfiltering.py and designfile.txt to current folder

USAGE: python qualityfiltering.py -d designfile.txt -f fastq.gz

1) It works with fastq.gz or fastq files
2) Changed the design format
3) Submits one job per sample
4) Checks if all the files (with user specified format) exists in the design file. If not the program exits without submitting the job
5) Time and memory options
6) Log file keeps record of sample to job id mapping
7) If qualityfiltered folder already exists, the program exits
'''

import sys,os
import shutil
import glob
from subprocess import call
import argparse as agp


def qualityfilteringPBS(designfile, fileformat):

        renamedict= dict()
        with open(designfile) as design:
                next(design)
                for val_d in design:
                        renamedict[val_d.split('\t')[0].strip()]= val_d.split('\t')[1].strip()
			renamedict[val_d.split('\t')[2].strip()]= val_d.split('\t')[3].strip()
        
	logfile= open('log_qualfilter.txt','wb')

        datadir= os.getcwd()
        print('Current working directory= ',datadir)
	for myfile in glob.glob('*.'+fileformat):
        	sbatch= open(datadir+'/sbatch_qualitytrim.sh','w')
        	sbatch.write('#!/bin/sh\n')
        	sbatch.write('#\n')
        	sbatch.write('#SBATCH --verbose\n')
        	sbatch.write('#SBATCH --job-name=qualityfilter\n')
        	sbatch.write('#SBATCH --output=qualfilter_%j.out\n')
        	sbatch.write('#SBATCH --error=qualfilter_%j.err\n')
        	sbatch.write('#SBATCH --time=01:50:00\n')
        	sbatch.write('#SBATCH --nodes=1\n')
        	sbatch.write('#SBATCH --cpus-per-task=1\n')
        	sbatch.write('#SBATCH --mem=10GB\n')
        	sbatch.write('\n')
        	sbatch.write('module load fastx_toolkit/intel/0.0.14\n')
        	sbatch.write('cd '+datadir+'\n')
        	sbatch.write('mkdir qualityfiltered\n')
                print('myfile= ',myfile)
                cmd_rename= ''.join(['mv ',myfile.strip(),' ',renamedict[myfile].strip()+'.'+fileformat])
                sbatch.write(cmd_rename+'\n')
                if '.gz' in fileformat: # check if file is zipped, unzip if zipped
                        sbatch.write('gunzip '+(renamedict[myfile].strip()+'.'+fileformat)+'\n')
                        myfile_unzip= (renamedict[myfile]+'.'+fileformat)[:-3]
                else:
                        myfile_unzip= renamedict[myfile]+'.'+fileformat
                trimclipout= 'qualityfiltered/'+renamedict[myfile].strip().split('.')[0]+'.clip.fq'
                mycmd1= ''.join(['fastx_clipper -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -l 20 -v -i ',myfile_unzip.strip(),' -o ',trimclipout,'\n'])
                sbatch.write(mycmd1)
                sbatch.write('\n')
        	sbatch.close()
        	shutil.copyfileobj(open("sbatch_qualitytrim.sh",'rb'), logfile)
        	call(["sbatch", "sbatch_qualitytrim.sh"])

if __name__=='__main__':
        parser= agp.ArgumentParser()
        parser.add_argument('-d',help='Provide design file for your data')
        parser.add_argument('-f',help='Provide your file format')
        args = parser.parse_args()
        qualityfilteringPBS(args.d,args.f)
		
