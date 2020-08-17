#!/usr/bin/env python3
import os
import pickle
import subprocess


def main(input_data):
    
    input_file = input_data[0]
    output_dir = input_data[1]
    ptf_path = input_data[2]
    translation_table = input_data[3]
    mode = input_data[4]

    if ptf_path != '':
        proc = subprocess.Popen(['prodigal', '-i', input_file, '-c',
                                 '-m', '-g', str(translation_table), '-p',
                                 mode, '-f', 'sco', '-q', '-t', ptf_path],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    elif ptf_path == '':
        proc = subprocess.Popen(['prodigal', '-i', input_file, '-c',
                                 '-m', '-g', str(translation_table), '-p',
                                 mode, '-f', 'sco', '-q'],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)

    # Read the stdout from Prodigal
    stdout = proc.stdout.readlines()
    stderr = proc.stderr.readlines()

    file_basename = os.path.basename(input_file).split('.')[0]

    if len(stderr) > 0:
        stderr = [line.decode('utf-8').strip() for line in stderr]
        stderr = [line for line in stderr if line != '']
        error = ' '.join(stderr)
        return [file_basename, error]

    # Parse output
    lines = [line.decode('utf-8').strip() for line in stdout]

    # determine contigs headers indexes
    contigs_headers = [l for l in lines if 'seqhdr' in l]
    contigs_ids = [l.split('"')[1].split()[0] for l in contigs_headers]
    contigs_idx = [lines.index(l) for l in contigs_headers] + [len(lines)]

    # get CDS' positions for each contig
    contigs_pos = {contigs_ids[i]: lines[contigs_idx[i]+1:contigs_idx[i+1]]
                   for i in range(len(contigs_ids))}

    # exclude contigs without coding sequences
    contigs_pos = {k: v[1:] for k, v in contigs_pos.items() if len(v) > 1}

    strand_trans = {'+': 1, '-': 0}

    # split and convert list elements
    contigs_pos = {k: [p.split('_')[1:] for p in v]
                   for k, v in contigs_pos.items()}
    contigs_pos = {k: [[int(p[0])-1, int(p[1]), strand_trans[p[2]]]
                   for p in v] for k, v in contigs_pos.items()}

    total_contigs = {k: len(v) for k, v in contigs_pos.items()}
    total_genome = sum(total_contigs.values())

    if total_genome > 0:
        # save positions in file
        filepath = os.path.join(output_dir, file_basename + '_ORF.txt')
        with open(filepath, 'wb') as f:
            pickle.dump(contigs_pos, f)

    status = [file_basename, total_genome]

    #print("done prodigal run on: " + str(os.path.basename(input_file)))

    return status


if __name__ == "__main__":

    main()
