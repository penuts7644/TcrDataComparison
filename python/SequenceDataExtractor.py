"""Extract CDR3 sequences with given V and J gene choices as well as full length
VDJ sequences from given CSV/TSV using the given IMGT genomic templates (species
and receptor type dependent).

"""

import argparse
import glob
import multiprocessing
import os
import sys

from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas
import numpy


def _parse_commandline():
    """Parses the commandline arguments given by user.

    Returns
    -------
    ArgumentParser
        Object containing parsed commandline arguments.

    """

    # Parse the commandline arguments given by user.
    parser = argparse.ArgumentParser(
        description="Extract CDR3 sequences with given V and J gene choices " \
                    "as well as full length VDJ sequences from given CSV/TSV " \
                    "using the given IMGT genomic templates (species and " \
                    "receptor type dependent)")
    parser.add_argument("files", metavar="F", type=str,
                        help="CSV or TSV file with sequence data or diretory with files")
    parser.add_argument("vgene", metavar="V", type=str,
                        help="V gene IMGT fasta file")
    parser.add_argument("jgene", metavar="J", type=str,
                        help="J gene IMGT fasts file")
    parser.add_argument("separator", metavar="S", type=str,
                        help="Separator charactor for input file")
    parser.add_argument("--num-threads", type=int, nargs="?", default=1,
                        help="Number of threads in integer (default: %(default)s)")
    parser.add_argument("--random-seqs", type=int, nargs="?",
                        help="Number of random sequences to use from each file " \
                             "in integer (default: all sequences)")
    return parser.parse_args()


def _read_fasta_as_dataframe(file):
    """Creates a pandas.DataFrame from the FASTA file.

    The dataframe contains header name and sequence columns containing the
    corresponding FASTA data.

    Parameters
    ----------
    file : string
        Location of the FASTA file to be read in.

    """
    # Create a dataframe and read in the fasta file.
    fasta_df = pandas.DataFrame(columns=['family', 'gene', 'nt_sequence'])
    with open(file, 'r') as fasta_file:
        for title, sequence in SimpleFastaParser(fasta_file):

            # Split up the gene identifier into seperate columns.
            full_id = title.split('|')[1].split('*')
            full_id[0] = full_id[0].split('-')
            if len(full_id[0]) == 2:
                family = full_id[0][0]
                gene = full_id[0][1]
            else:
                family = full_id[0][0]
                gene = numpy.nan

            # Add row to DataFrame.
            fasta_df = fasta_df.append({
                'family': family,
                'gene': gene,
                'nt_sequence': ''.join(sequence.upper().split('.')),
            }, ignore_index=True)
    return fasta_df


def _find_longest_substring(full, partial):
    """Counts the combined occurrences of each Kmer for pandas.DataFrame rows.

    Parameters
    ----------
    full : str
        A full length sequence string.
    partial : str
        A partial length sequence string to compare to the ful length one.

    Returns
    -------
    str
        The longest substring from the compared input strings.

    """
    t = [[0] * (1 + len(partial)) for i in range(1 + len(full))]
    l, xl = 0, 0
    for x in range(1, 1 + len(full)):
        for y in range(1, 1 + len(partial)):
            if full[x - 1] == partial[y - 1]:
                t[x][y] = t[x - 1][y - 1] + 1
                if t[x][y] > l:
                    l = t[x][y]
                    xl  = x
            else:
                t[x][y] = 0
    return full[xl - l: xl]


def reassemble_data(args):
    """Reassembles the dataframe data by striping the CDR3 and creating the VDJ.

    Parameters
    ----------
    args : list
        A collection of arguments containing dataframe, V and J gene datframes
        as well as the correct column names to use.

    Returns
    -------
    DataFrames
        Four pandas datframes containing the reassembled data, full length
        productive VDJ sequences, full length unproductive VDJ sequences and
        one with all full length VDJ sequences.

    Notes
    -----
    The gene name is separated into three sections: family, gene and allele. By
    default, family and gene are required while for the allele '*01' is used.

    """
    # Setup the initial dataframe.
    df, kwargs = args
    reassembled_df = pandas.DataFrame(columns=[
        'seq_id', 'nt_sequence', 'aa_sequence', 'gene_choice_v', 'gene_choice_j'
    ])
    full_length_prod_df = pandas.DataFrame(columns=['seq_id', 'nt_sequence'])
    full_length_unprod_df = pandas.DataFrame(columns=['seq_id', 'nt_sequence'])
    full_length_df = pandas.DataFrame(columns=['seq_id', 'nt_sequence'])

    # Iterate over the rows with index value.
    for i, row in df.iterrows():

        gene_choice_v = numpy.nan
        imgt_v_gene = pandas.DataFrame()
        if (isinstance(row[kwargs['col_names']['v_family']], str)
                and isinstance(row[kwargs['col_names']['v_gene']], str)):

            # Pre-process the V gene.
            v_family = row[kwargs['col_names']['v_family']].replace('TCR', 'TR').replace('V0', 'V')
            v_gene = row[kwargs['col_names']['v_gene']].split('-')[1].split('/')[0].lstrip('0')
            imgt_v_gene = kwargs['v_genes'][
                (kwargs['v_genes']['family'] == v_family)
                & ((kwargs['v_genes']['gene'] == v_gene)
                   | ((kwargs['v_genes']['gene'].isna()) & (v_gene == '1')))].head(1)

            # Assemble the output dataframe gene choices for the V gene.
            if not imgt_v_gene.empty:
                if not isinstance(imgt_v_gene['gene'].values[0], str) and v_gene == '1':
                    gene_choice_v = imgt_v_gene['family'].values[0] + '*01'
                else:
                    gene_choice_v = imgt_v_gene['family'].values[0] + '-' \
                                    + imgt_v_gene['gene'].values[0] + '*01'

        gene_choice_j = numpy.nan
        imgt_j_gene = pandas.DataFrame()
        if (isinstance(row[kwargs['col_names']['j_family']], str)
                and isinstance(row[kwargs['col_names']['j_gene']], str)):

            # Pre-process the J gene.
            j_family = row[kwargs['col_names']['j_family']].replace('TCR', 'TR').replace('J0', 'J')
            j_gene = row[kwargs['col_names']['j_gene']].split('-')[1].split('/')[0].lstrip('0')
            imgt_j_gene = kwargs['j_genes'][
                (kwargs['j_genes']['family'] == j_family)
                & ((kwargs['j_genes']['gene'] == j_gene)
                   | ((kwargs['j_genes']['gene'].isna()) & (j_gene == '1')))].head(1)

            # Assemble the output dataframe gene choices for the J gene.
            if not imgt_j_gene.empty:
                if not isinstance(imgt_j_gene['gene'].values[0], str) and j_gene == '1':
                    gene_choice_j = imgt_j_gene['family'].values[0] + '*01'
                else:
                    gene_choice_j = imgt_j_gene['family'].values[0] + '-' \
                                    + imgt_j_gene['gene'].values[0] + '*01'

        # Create the trimmed NT sequence (removing primers).
        trimmed_nt_seq = row[kwargs['col_names']['nt_seq']][
            (81 - int(row[kwargs['col_names']['cdr3_len']])): 81]

        # Add data row of reassembled data to the dataframe.
        reassembled_df = reassembled_df.append({
            'seq_id': i,
            'nt_sequence': trimmed_nt_seq,
            'aa_sequence': row[kwargs['col_names']['aa_seq']],
            'gene_choice_v': gene_choice_v,
            'gene_choice_j': gene_choice_j,
        }, ignore_index=True)

        # Create the VDJ full length sequence
        if (not imgt_v_gene.empty
                and not imgt_j_gene.empty):
            vd_segment = _find_longest_substring(
                imgt_v_gene['nt_sequence'].values[0], trimmed_nt_seq)
            dj_segment = _find_longest_substring(
                imgt_j_gene['nt_sequence'].values[0], trimmed_nt_seq)
            split_v = imgt_v_gene['nt_sequence'].values[0].rsplit(vd_segment, 1)
            split_j = imgt_j_gene['nt_sequence'].values[0].split(dj_segment, 1)
            if (len(split_v[1]) >= len(split_v[0])) or (len(split_j[0]) >= len(split_j[1])):
                continue
            vdj_sequence = split_v[0] + trimmed_nt_seq + split_j[1]
        else:
            continue

        # Add data row of full length data to the dataframe for productive and
        # unproductive sequences.
        if row[kwargs['col_names']['type']].lower() == 'in':
            full_length_prod_df = full_length_prod_df.append({
                'seq_id': i,
                'nt_sequence': vdj_sequence
            }, ignore_index=True)
        elif (row[kwargs['col_names']['type']].lower() == 'out'
              or row[kwargs['col_names']['type']].lower() == 'stop'):
            full_length_unprod_df = full_length_unprod_df.append({
                'seq_id': i,
                'nt_sequence': vdj_sequence
            }, ignore_index=True)
        full_length_df = full_length_df.append({
            'seq_id': i,
            'nt_sequence': vdj_sequence
        }, ignore_index=True)
    return reassembled_df, full_length_prod_df, full_length_unprod_df, full_length_df


def multiprocess_dataframe(df, func, num_workers, **kwargs):
    """Applies multiprocessing on a pandas dataframe using given function.

    Parameters
    ----------
    df : DataFrame
        Pandas dateframe to be split for multiple workers.
    func : Object
        A function object that the workers should apply.
    num_workers : int
        Integer specifying the number of workers (threads) to create.
    **kwargs : dict
        The remaining arguments to be given to the input function.

    Returns
    -------
    list
        Contains the results from each of the workers.

    """
    # Check out available worker count and adjust accordingly.
    if len(df) < num_workers:
        num_workers = len(df)

    # Divide the array into chucks for the workers.
    pool = multiprocessing.Pool(processes=num_workers)
    result = pool.map(func, [(d, kwargs)
                             for d in numpy.array_split(df, num_workers)])
    pool.close()
    return result


def main():
    """Main function when script ran through command-line."""

    # Set some default parameters for the input file.
    column_names = {
        "type": "frame_type",
        "nt_seq": "rearrangement",
        "aa_seq": "amino_acid",
        "cdr3_len": "cdr3_length",
        "v_family": "v_family",
        "v_gene": "v_gene",
        "j_family": "j_family",
        "j_gene": "j_gene",
    }

    # Parse arguments and read in given files.
    args = _parse_commandline()

    # Locate the files to process.
    files = None
    if os.path.isdir(args.files):
        if r'\t' in args.separator:
            files = [f for f in glob.glob(os.path.join(args.files, '*.tsv'))]
        else:
            files = [f for f in glob.glob(os.path.join(args.files, '*.csv'))]
    elif os.path.isfile(args.files):
        files = [args.files]
    if not files:
        print('No files where found in the given directory')
        sys.exit()

    vgene_df = _read_fasta_as_dataframe(args.vgene)
    jgene_df = _read_fasta_as_dataframe(args.jgene)

    reassembled_df = pandas.DataFrame()
    full_length_prod_df = pandas.DataFrame()
    full_length_unprod_df = pandas.DataFrame()
    full_length_df = pandas.DataFrame()

    for file in files:
        data_df = pandas.read_csv(
            file, sep=args.separator, header=0, comment='#', dtype=str,
            na_values=['na', 'unknown', 'unresolved', 'no data'], engine='python',
            usecols=list(column_names.values())
        )

        # Take a random subsample of sequences in the file.
        if args.random_seqs is not None:
            if args.random_seqs > 0 and len(data_df) >= args.random_seqs:
                data_df = data_df.sample(n=args.random_seqs, random_state=1)
            else:
                print('Number of random sequences should be higher 0 and smaller ' \
                      'than total number of rows')
                sys.exit()

        # Process the dataframe in various threads and append results.
        results = multiprocess_dataframe(
            df=data_df, func=reassemble_data, num_workers=args.num_threads,
            v_genes=vgene_df, j_genes=jgene_df, col_names=column_names)
        for processed in results:
            processed[0].insert(0, 'file_name', os.path.splitext(os.path.basename(file))[0])
            processed[1].insert(0, 'file_name', os.path.splitext(os.path.basename(file))[0])
            processed[2].insert(0, 'file_name', os.path.splitext(os.path.basename(file))[0])
            processed[3].insert(0, 'file_name', os.path.splitext(os.path.basename(file))[0])
            reassembled_df = reassembled_df.append(processed[0], ignore_index=True)
            full_length_prod_df = full_length_prod_df.append(processed[1], ignore_index=True)
            full_length_unprod_df = full_length_unprod_df.append(processed[2], ignore_index=True)
            full_length_df = full_length_df.append(processed[3], ignore_index=True)

    # Writes the new dataframe to a CSV file.
    output_filename_base = 'sequence_data_extractor'
    output_filename_1 = os.path.join(os.getcwd(), output_filename_base + "_CDR3.tsv")
    # reassembled_df.reset_index(inplace=True)
    pandas.DataFrame.to_csv(reassembled_df, path_or_buf=output_filename_1,
                            sep="\t", na_rep="NA", index=True, index_label='seq_index')
    print("Written '{}' file".format(output_filename_1))
    output_filename_2 = os.path.join(os.getcwd(), output_filename_base + "_productive.tsv")
    # full_length_prod_df.reset_index(inplace=True)
    pandas.DataFrame.to_csv(full_length_prod_df, path_or_buf=output_filename_2,
                            sep="\t", na_rep="NA", index=True, index_label='seq_index')
    print("Written '{}' file".format(output_filename_2))
    output_filename_3 = os.path.join(os.getcwd(), output_filename_base + "_unproductive.tsv")
    # full_length_unprod_df.reset_index(inplace=True)
    pandas.DataFrame.to_csv(full_length_unprod_df, path_or_buf=output_filename_3,
                            sep="\t", na_rep="NA", index=True, index_label='seq_index')
    print("Written '{}' file".format(output_filename_3))
    output_filename_4 = os.path.join(os.getcwd(), output_filename_base + "_all.tsv")
    # full_length_df.reset_index(inplace=True)
    pandas.DataFrame.to_csv(full_length_df, path_or_buf=output_filename_4,
                            sep="\t", na_rep="NA", index=True, index_label='seq_index')
    print("Written '{}' file".format(output_filename_4))


if __name__ == "__main__":
    main()
