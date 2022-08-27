#!/usr/bin/python

''' 

MACS2 by default removes the duplicates from your file before peak calling. 
However, it provides an option if you want to consider duplicate reads. Provide a
parameter -f withdup to the python script, if you want to run MACS2 with duplicate reads.
For ChIP-seq and DamID-seq we removed the duplicates. 

Copy runmacs2.py and designfile.txt to current folder.

USAGE: python runmacs2.py -d designfile.txt -f nodup
 
'''

import sys,os
import shutil
import glob
from subprocess import call
import argparse as agp

def domapping(designfile,duplicates):

        chip_control_dict= dict()
        with open(designfile) as design:
                next(design)
                for val_d in design:
                        chip_control_dict[val_d.split('\t')[1].split('.')[0].strip()]= val_d.split('\t')[3].split('.')[0].strip()
	#print('chip_control_dict= ',chip_control_dict)
	
	logfile= open('log_peakcall.txt','wb')	
	
        datadir= os.getcwd()
        print 'Current working directory= ',datadir
        sbatch= open(datadir+'/sbatch_macs2.sh','w')
        sbatch.write('#!/bin/sh\n')
        sbatch.write('#\n')
        sbatch.write('#SBATCH --verbose\n')
        sbatch.write('#SBATCH --job-name=run_macs2\n')
        sbatch.write('#SBATCH --output=macs2_%j.out\n')
        sbatch.write('#SBATCH --error=macs2_%j.err\n')
        sbatch.write('#SBATCH --time=02:00:00\n')
        sbatch.write('#SBATCH --nodes=1\n')
        sbatch.write('#SBATCH --cpus-per-task=4\n')
        sbatch.write('#SBATCH --mem=20GB\n')
        sbatch.write('\n')
        sbatch.write('cd '+datadir+'\n')
        sbatch.write('module load macs2/intel/2.1.1\n')
        sbatch.write('# Run macs2\n')
	for val_dict in chip_control_dict:
		if duplicates.strip()=='nodup' or duplicates.strip()==None:
                	macs2cmd= ''.join(['macs2 callpeak -t ',(val_dict+'.nodup.bam'),' -c ',(chip_control_dict[val_dict]+'.nodup.bam'),' -f BAM -g 1.2e8 -n ',(val_dict+'.macs2 -q 0.05'),'\n'])
                if duplicates.strip()=='withdup':
			macs2cmd= ''.join(['macs2 callpeak -t ',(val_dict+'.sorted.bam'),' -c ',(chip_control_dict[val_dict]+'.sorted.bam'),' -f BAM -g 1.2e8 --keep-dup -n ',(val_dict+'.macs2 -q 0.05'),'\n'])
		sbatch.write(macs2cmd)
	sbatch.write('\n')
	sbatch.close()
	shutil.copyfileobj(open("sbatch_macs2.sh",'rb'), logfile)
	call(["sbatch", "sbatch_macs2.sh"])

if __name__=='__main__':
        parser= agp.ArgumentParser()
        parser.add_argument('-d',help='Provide design file for your data')
        #parser.add_argument('-p',help='Provide path for scripts')
	parser.add_argument('-f',help='Options:nodup or withdup',default='nodup')
        args = parser.parse_args()
        domapping(args.d,args.f)
