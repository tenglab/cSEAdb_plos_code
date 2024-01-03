#!/usr/bin/env python
"""convert bowtie output sam to sort bam"""

import argparse
import os
import sys

# --------------------------------------------------
def get_args():
	"""get args"""
	parser = argparse.ArgumentParser(description='count bowtie sam')
	parser.add_argument('-i', '--infile', help='bowtie sam',
						type=str, metavar='FILE', required=True)
	parser.add_argument('-od', '--out_dir', help='output dir',
						type=str, metavar='PATH', required=True)	
	parser.add_argument('-o', '--outfile', help='output file',
						type=str, metavar='FILE', required=True)
	return parser.parse_args()

# --------------------------------------------------
def main():
	"""main"""
	args = get_args()
	infile = args.infile
	out_dir = args.out_dir
	outfile = args.outfile

	in_dir = os.path.dirname(os.path.abspath(infile))
	temp_bam = out_dir+'/temp_'+outfile+'.bam'	
	temp_sort_bam = out_dir+'/temp_'+outfile+'_sorted.bam'
	final_bam = out_dir+'/'+outfile+'_rmdup_sort.bam'

	os.chdir(in_dir)
	os.system('samtools view -bS '+infile+'> '+temp_bam)
	os.system('samtools sort '+temp_bam+' -o '+temp_sort_bam)
	os.system('samtools rmdup '+temp_sort_bam+' '+final_bam)
	os.system('samtools index '+final_bam)
	os.system('rm -rf '+temp_bam)
	os.system('rm -rf '+temp_sort_bam)

# --------------------------------------------------
if __name__ == '__main__':
	main()
