#!/usr/bin/env python3


import csv
import argparse


def main(input_file, genes_list, output_file, inverse):

    if inverse:
        FilesToRemove = ['File', 'FILE', 'file']
    else:
        FilesToRemove = []

    with open(genes_list) as f:
        i = 0
        for File in f:
            File = File.rstrip('\n')
            File = File.rstrip('\r')
            File = (File.split('\t'))[0]
            FilesToRemove.append(File)
            i += 1

    print('\nProvided list has {0} genes.'.format(i-1))

    with open(input_file, 'r') as tsvin, open(output_file + ".tsv", "w") as csvout:
        tsvin = csv.reader(tsvin, delimiter='\t')

        listindextoremove = []
        for firstline in tsvin:
            for gene in firstline:
                if gene in FilesToRemove and not inverse:
                    listindextoremove.append(firstline.index(gene))
                elif gene not in FilesToRemove and inverse:
                    listindextoremove.append(firstline.index(gene))

            print('Removing {0} genes...'.format(len(listindextoremove)), end='')

            for elem in reversed(listindextoremove):
                del firstline[elem]
            csvout.write(('\t'.join(firstline)) + "\n")
            break

        for line in tsvin:
            for elem in reversed(listindextoremove):
                del line[elem]
            csvout.write(('\t'.join(line)) + "\n")

        print('done.')


if __name__ == "__main__":

    main()
