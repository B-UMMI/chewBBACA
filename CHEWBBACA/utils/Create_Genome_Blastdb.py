#!/usr/bin/env python3


import os
import sys


def main(makeblastdb_path, questionDB, directory, genomeFile, nucleotide=False):

	name = directory + "/" + genomeFile + "_db"

	cmd_template = '{0} -in {1} -out {2} -dbtype {3} -logfile {2}_blast.log'

	if nucleotide is True:
		os.system(cmd_template.format(makeblastdb_path, questionDB, name, 'nucl'))
	else:
		os.system(cmd_template.format(makeblastdb_path, questionDB, name, 'prot'))

	return True


if __name__ == "__main__":

	main()
