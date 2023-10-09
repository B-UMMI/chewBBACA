#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This module generates an interactive HTML report for the allele calling
results. The report provides summary statistics to evaluate results per
sample and per locus (with the possibility to provide a TSV file with
loci annotations to include on a table). The report includes components
to display a heatmap representing the loci presence-absence matrix, a
heatmap representing the distance matrix based on allelic differences
and a Neighbor-Joining (NJ) tree based on the MSA of the core genome loci.

Code documentation
------------------
"""


import os
import json
import pandas as pd

try:
    from utils import (
        constants as ct,
        file_operations as fo,
        fasta_operations as fao,
        iterables_manipulation as im,
        multiprocessing_operations as mo,
        distance_matrix as dm,
        mafft_wrapper as mw,
        fasttree_wrapper as fw)
    from ExtractCgMLST import determine_cgmlst
except ModuleNotFoundError:
    from CHEWBBACA.utils import (
        constants as ct,
        file_operations as fo,
        fasta_operations as fao,
        iterables_manipulation as im,
        multiprocessing_operations as mo,
        distance_matrix as dm,
        mafft_wrapper as mw,
        fasttree_wrapper as fw)
    from CHEWBBACA.ExtractCgMLST import determine_cgmlst


def compute_sample_stats(sample_ids, total_loci, coordinates_file,
                         sample_counts):
    """Compute sample statistics to include in the report.

    Parameters
    ----------
    sample_ids : list
        List of sample unique identifiers.
    total_loci : list
        List of loci unique identifiers.
    coordinates_file : str
        Path to the TSV file with CDS coordinates data created by
        the AlleleCall module.
    sample_counts : pandas.core.frame.DataFrame
        Dataframe with the data from the 'results_statistics.tsv' file
        created by the AlleleCall module.

    Returns
    -------
    sample_stats : dict
        Dictionary with sample unique identifiers as keys and a list
        with the computed statistics as values (number of contigs,
        number of CDSs, list of contig identifiers, proportion of
        CDSs classified, number of identified loci, proportion of
        identified loci, valid classifications, invalid classifications
        and the count for each classification type).
    """
    sample_stats = {sid: [0, 0, []] for sid in sample_ids}
    # Get number of contigs, CDSs and list of contig ids per sample
    with open(coordinates_file, 'r') as infile:
        # Skip header line
        header = infile.__next__()
        for line in infile:
            current_line = line.strip().split('\t')
            if current_line[1] not in sample_stats[current_line[0]][2]:
                sample_stats[current_line[0]][2].append(current_line[1])
                sample_stats[current_line[0]][0] += 1
            sample_stats[current_line[0]][1] += 1

    # Get class counts per sample
    for i, sample in enumerate(sample_stats):
        sample_line = sample_counts.iloc[i].tolist()
        # Compute stats
        # EXC + INF
        valid = sum(sample_line[1:3])
        # PLOT3 + PLOT5 + LOTSC + NIPH + NIPHEM + ALM + ASM + PAMALNF"
        invalid = sum(sample_line[3:11])
        total_cds = sample_stats[sample][1]
        # Only counts 1 CDS per NIPH/NIPHEM
        identified_loci = valid + invalid
        classified_cds_proportion = round(identified_loci/total_cds, 3)
        identified_loci_proportion = round(identified_loci/total_loci, 3)
        # Add computed values to sample stats
        sample_stats[sample].append(float(classified_cds_proportion))
        sample_stats[sample].append(int(identified_loci))
        sample_stats[sample].append(float(identified_loci_proportion))
        sample_stats[sample].append(int(valid))
        sample_stats[sample].append(int(invalid))
        # Add class counts
        sample_stats[sample].extend(list(map(int, sample_line[1:])))

    return sample_stats


def compute_loci_stats(total_samples, loci_counts):
    """Compute loci statistics to include in the report.

    Parameters
    ----------
    total_samples : list
        Total number of samples.
    loci_counts: pandas.core.frame.DataFrame
        Dataframe with the data from the 'loci_summary_stats.tsv' file
        created by the AlleleCall module.

    Returns
    -------
    loci_stats : list
        List with one sublist per locus. Each sublist includes a
        locus unique identifier, the total number of CDSs classified
        for that locus, the number of valid classifications, the
        number of invalid classifications, the proportion of samples
        the locus was found in and the count for each classification type.
    """
    loci_stats = []
    for i in range(len(loci_counts)):
        locus_line = loci_counts.iloc[i]
        locus_id = locus_line['Locus']
        total_cds = int(locus_line['Total_CDS'])
        valid = sum(locus_line[1:3])
        invalid = sum(locus_line[3:11])
        sample_proportion = round(valid/total_samples, 3)
        loci_stats.append([locus_id, int(total_cds), int(valid),
                           int(invalid), float(sample_proportion),
                           *list(map(int, locus_line[1:12]))])

    return loci_stats


def profile_column_to_fasta(locus, schema_directory, allelic_profiles,
                            output_directory):
    """Create a FASTA file with the locus alleles identified in the dataset.

    Parameters
    ----------
    locus : str
        The locus identifier.
    schema_directory : str
        Path to the schema directory.
    allelic_profiles : str
        Path to the 'results_alleles.tsv' file create by the AlleleCall
        module.
    output_directory : str
        Path to the output directory.

    Returns
    -------
    fasta_file : str
        Path to the FASTA file that contains the allele sequences
        identified for each sample.
    """
    # Read locus column
    df = pd.read_csv(allelic_profiles, usecols=['FILE', locus],
                     delimiter='\t', dtype=str)
    locus_column = df[locus].tolist()
    sample_ids = df['FILE'].tolist()

    # Read locus FASTA file
    locus_fasta = fo.join_paths(schema_directory, [f'{locus}.fasta'])
    alleles = fao.import_sequences(locus_fasta)
    alleles = {k.split('_')[-1]: v for k, v in alleles.items()}

    sequences = []
    for i, allele in enumerate(locus_column):
        clean_class = allele.split('INF-')[-1]
        if clean_class in alleles:
            current_sample = sample_ids[i]
            seq = alleles[clean_class]
            record = fao.fasta_str_record(ct.FASTA_RECORD_TEMPLATE,
                                          [f'{locus}_{current_sample}', seq])
            sequences.append(record)

    fasta_file = fo.join_paths(output_directory, [f'{locus}.fasta'])
    fo.write_lines(sequences, fasta_file)

    return fasta_file


def concatenate_loci_alignments(sample, loci, fasta_index, output_directory):
    """Concatenate the aligned sequences for a sample.

    Parameters
    ----------
    sample : str
        Sample identifier.
    loci : list
        Loci identifiers.
    fasta_index : Bio.File._IndexedSeqFileDict
        Indexed FASTA file to get sequences from.
    output_directory : str
        Path to the output directory.
    """
    alignment = ''
    for locus in loci:
        # Sequence headers include locus and sample IDs joined by '_'
        seqid = f'{locus}_{sample}'
        # Get aligned sequence from index
        try:
            alignment += str(fasta_index[seqid].seq)
        except Exception as e:
            print(f'Could not get {sample} allele for locus {locus}.')
    # Save alignment for sample
    alignment_outfile = fo.join_paths(output_directory,
                                      [f'{sample}_cgMLST_alignment.fasta'])
    alignment_record = fao.fasta_str_record(ct.FASTA_RECORD_TEMPLATE,
                                            [sample, alignment])
    fo.write_lines([alignment_record], alignment_outfile)

    return alignment_outfile


def main(input_files, schema_directory, output_directory, annotations,
         cpu_cores, light, no_pa, no_dm, no_tree, cg_alignment):

    # Create temp directory
    temp_directory = fo.join_paths(output_directory, ['temp'])
    fo.create_directory(temp_directory)

    # Read "results_statistics.tsv" to get class counts per sample
    sample_statistics_file = fo.join_paths(input_files,
                                           [ct.RESULTS_STATISTICS_BASENAME])
    sample_counts = pd.read_csv(sample_statistics_file, delimiter='\t')

    # Sort based on decreasing number of EXC
    sample_counts = sample_counts.sort_values(by=['EXC'], ascending=False)

    # Get sample identifiers
    sample_ids = sample_counts['FILE'].tolist()
    # Get total number of samples
    total_samples = len(sample_ids)
    print(f'Number of samples: {total_samples}')

    # Data for stacked Bar chart with class counts per sample
    sample_data = []
    for c in ct.ALLELECALL_CLASSIFICATIONS:
        sample_data.append(sample_counts[c].tolist())

    # Read "loci_summary_stats.tsv" to get class counts per locus
    loci_statistics_file = fo.join_paths(input_files,
                                         [ct.LOCI_STATS_BASENAME])
    loci_counts = pd.read_csv(loci_statistics_file, delimiter='\t')
    # Sort based on decreasing number of EXC
    loci_counts = loci_counts.sort_values(by=['EXC'], ascending=False)

    # Get loci identifiers
    loci_ids = loci_counts['Locus'].tolist()
    # Get total number of loci
    total_loci = len(loci_ids)
    print(f'Number of loci: {total_loci}')

    # Data for stacked Bar chart with class counts per locus
    loci_data = []
    for c in ct.ALLELECALL_CLASSIFICATIONS:
        loci_data.append(loci_counts[c].tolist())

    # Compute data for Sample Stats table
    # Need to read file with CDS coordinates to get number of contigs and
    # CDSs per sample (future development will use coordinates for synteny
    # analysis)
    print('Computing sample statistics...', end='')
    cds_coordinates_file = fo.join_paths(input_files, [ct.CDS_COORDINATES_BASENAME])
    sample_stats = compute_sample_stats(sample_ids, total_loci,
                                        cds_coordinates_file, sample_counts)

    total_cds = sum([v[1] for k, v in sample_stats.items()])

    sample_stats_columns = ct.SAMPLE_STATS_COLUMNS + \
        sample_counts.columns.tolist()[1:]

    # Create list to store in HTML
    sample_stats_rows = []
    for k, v in sample_stats.items():
        sample_stats_rows.append([k, v[0], v[1], *v[3:]])
    print('done.')

    # Compute data for the Loci Stats table
    print('Computing loci statistics...', end='')
    loci_stats = compute_loci_stats(total_samples, loci_counts)
    # Sort loci stats based on decreasing sample presence
    loci_stats = sorted(loci_stats, key=lambda x: x[4], reverse=True)

    loci_stats_columns = ct.LOCI_STATS_COLUMNS + \
        loci_counts.columns.tolist()[1:-1]
    print('done.')

    # Read loci annotations from TSV file
    annotation_values = [[], []]
    if annotations is not None:
        annotation_lines = fo.read_tabular(annotations)
        # Add file columns
        annotation_values[0].extend(annotation_lines[0])
        # Add sublists with lines
        # Only keep lines for loci in the schema
        loci_annotations = [line for line in annotation_lines[1:]
                            if line[0] in loci_ids]
        annotation_values[1].extend(loci_annotations)

    print(f'Provided annotations for {len(annotation_values[1])} '
          'loci in the schema.')

    # Data for Summary Table
    # Count total number of classified CDSs per class
    loci_sums = loci_counts[loci_counts.columns[1:]].sum(axis=0)
    loci_sums = loci_sums.tolist()
    summary_columns = ct.SUMMARY_STATS_COLUMNS
    summary_rows = [total_samples, total_loci, total_cds,
                    loci_sums[-1], *loci_sums[:-1]]

    pa_lines = []
    dm_lines = []
    phylo_data = {"phylo_data": []}
    if light is False:
        if False in [no_pa, no_dm, no_tree] or cg_alignment is True:
            # Define path to TSV file that contains allelic profiles
            allelic_profiles_file = fo.join_paths(input_files,
                                                  [ct.RESULTS_ALLELES_BASENAME])

            # Import matrix with allelic profiles
            print('Reading profile matrix...', end='')
            profiles_matrix = pd.read_csv(allelic_profiles_file,
                                          header=0, index_col=0,
                                          sep='\t', low_memory=False)
            print('done.')
            # Mask missing data
            print('Masking profile matrix...', end='')
            masked_profiles = profiles_matrix.apply(im.replace_chars)
            output_masked = os.path.join(output_directory, ct.MASKED_PROFILES_BASENAME)
            masked_profiles.to_csv(output_masked, sep='\t')
            print('done.')
            # Compute Presence-Absence matrix
            print('Computing Presence-Absence matrix...', end='')
            pa_matrix, pa_outfile = determine_cgmlst.presAbs(masked_profiles,
                                                             output_directory)
            print('done.')

            if no_pa is False:
                # Sort Presence-Absence matrix based on decreasing loci presence
                sorted_loci = [x[0] for x in loci_stats]
                pa_matrix = pa_matrix[sorted_loci]
                sorted_samples = pa_matrix.index.tolist()
                pa_data = [{"rows": pa_matrix.values.tolist()},
                           {"loci_ids": sorted_loci},
                           {"sample_ids": sorted_samples}]

            if no_dm is False or no_tree is False or cg_alignment is True:
                # Compute the cgMLST at 100%
                print('Determining cgMLST loci...')
                cgMLST_genes, _ = determine_cgmlst.compute_cgMLST(pa_matrix, sample_ids,
                                                                  1, len(sample_ids))
                cgMLST_genes = cgMLST_genes.tolist()
                print('\n', f'cgMLST is composed of {len(cgMLST_genes)} loci.')
                cgMLST_matrix = masked_profiles[cgMLST_genes]
                cgMLST_matrix_outfile = os.path.join(output_directory, ct.CGMLST_PROFILES_BASENAME)
                cgMLST_matrix.to_csv(cgMLST_matrix_outfile, sep='\t')

            if no_dm is False:
                # Compute distance matrix
                # Based on cgMLST profiles
                dm_file = dm.main(cgMLST_matrix_outfile, output_directory,
                                  cpu_cores, True, True)
                # Import distance matrix
                distance_m = pd.read_csv(dm_file[0], header=0, index_col=0,
                                         sep='\t', low_memory=False)
                dm_data = [{"rows": distance_m.values.tolist()},
                           {"sample_ids": distance_m.columns.tolist()}]

        # Only using the loci in the cgMLST
        # Might have to change if we need to work with all loci in the future
        if no_tree is False or cg_alignment is True:
            # Create FASTA files with alleles identified in samples
            # Create temporary directory to store FASTA files
            print('Creating FASTA files with identified alleles...')
            fasta_dir = fo.join_paths(temp_directory, ['fasta_files'])
            fo.create_directory(fasta_dir)

            inputs = im.divide_list_into_n_chunks(cgMLST_genes, len(cgMLST_genes))

            common_args = [schema_directory, allelic_profiles_file, fasta_dir]

            # Add common arguments to all sublists
            inputs = im.multiprocessing_inputs(inputs,
                                               common_args,
                                               profile_column_to_fasta)

            # Create FASTA files with identified alleles
            results = mo.map_async_parallelizer(inputs,
                                                mo.function_helper,
                                                cpu_cores,
                                                show_progress=True)

            # Translate FASTA files
            print('\nTranslating FASTA files...')
            translation_inputs = im.divide_list_into_n_chunks(results, len(results))
            common_args = [fasta_dir, ct.GENETIC_CODES_DEFAULT]
            # Add common arguments to all sublists
            inputs = im.multiprocessing_inputs(translation_inputs,
                                               common_args,
                                               fao.translate_fasta)
            results = mo.map_async_parallelizer(inputs,
                                                mo.function_helper,
                                                cpu_cores,
                                                show_progress=True)

            protein_files = [r[1] for r in results]

            # Align sequences with MAFFT
            print('\nDetermining the MSA for each locus...')
            alignment_dir = fo.join_paths(temp_directory, ['alignment_files'])
            fo.create_directory(alignment_dir)
            alignmnet_inputs = im.divide_list_into_n_chunks(protein_files,
                                                            len(protein_files))
            common_args = [alignment_dir]

            # Add common arguments to all sublists
            inputs = im.multiprocessing_inputs(alignmnet_inputs,
                                               common_args,
                                               mw.call_mafft)

            results = mo.map_async_parallelizer(inputs,
                                                mo.function_helper,
                                                cpu_cores,
                                                show_progress=True)
            mafft_files = {fo.file_basename(file, False).split('_protein_aligned')[0]: file
                           for file in results if file is not False}

            print('\nCreating file with the full cgMLST alignment...', end='')
            # Concatenate all alignment files and index with BioPython
            concat_aln = fo.join_paths(alignment_dir, ['cgMLST_concat.fasta'])
            fo.concatenate_files(mafft_files.values(), concat_aln)
            # Index file
            concat_index = fao.index_fasta(concat_aln)
            sample_alignment_files = []
            # Not dealing well with '*' in allele ids
            for sample in sample_ids:
                alignment_file = concatenate_loci_alignments(sample,
                                                             cgMLST_genes,
                                                             concat_index,
                                                             fasta_dir)
                sample_alignment_files.append(alignment_file)

            # Concatenate all cgMLST alignmnet records
            full_alignment = fo.join_paths(output_directory,
                                           [ct.CORE_MSA_BASENAME])
            fo.concatenate_files(sample_alignment_files, full_alignment)
            print('done.')

            if no_tree is False:
                print('Computing the NJ tree based on the core genome MSA...', end='')
                # Compute NJ tree with FastTree
                out_tree = fo.join_paths(alignment_dir, ['cgMLST.tree'])
                fw.call_fasttree(full_alignment, out_tree)
                phylo_data = fo.read_file(out_tree)
                phylo_data = {"phylo_data": phylo_data}
                print('done.')

    report_data = {"summaryData": [{"columns": summary_columns},
                                   {"rows": [summary_rows]}],
                   "sample_ids": sample_ids,
                   "sample_data": sample_data,
                   "loci_ids": loci_ids,
                   "loci_data": loci_data,
                   "sample_stats": [{"columns": sample_stats_columns},
                                    {"rows": sample_stats_rows}],
                   "loci_stats": [{"columns": loci_stats_columns},
                                  {"rows": loci_stats}],
                   "annotations": [{"columns": annotation_values[0]},
                                   {"rows": annotation_values[1]}],
                   "presence_absence": pa_data,
                   "distance_matrix": dm_data,
                   "cgMLST_tree": phylo_data,
                   "cgMLST_size": [{"cgMLST100": len(cgMLST_genes)}],
                   }

    # Write report HTML file
    report_html = ct.ALLELECALL_REPORT_HTML
    report_html = report_html.format(json.dumps(report_data))
    report_html_file = fo.join_paths(output_directory, [ct.ALLELECALL_REPORT_BASENAME])
    fo.write_to_file(report_html, report_html_file, 'w', '\n')

    # Copy JS bundle to output directory
    script_path = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(script_path)
    try:
        fo.copy_file(fo.join_paths(script_path,
                                   ['report_template_components',
                                    'src',
                                    'bundles',
                                    'AlleleCallEvaluator',
                                    'report_bundle.js']),
                     output_directory)
    except Exception as e:
        fo.copy_file(fo.join_paths(parent_dir,
                                   ['report_template_components',
                                    'src',
                                    'bundles',
                                    'AlleleCallEvaluator',
                                    'report_bundle.js']),
                     output_directory)

    # Delete all temporary files
    fo.delete_directory(temp_directory)

    print(f'\nResults available in {output_directory}.')
