#!/usr/bin/env python

"""
Requirements: Biopython

This script takes an input fastq file (compressed with gzip, i.e. suffix .gz, or uncompressed), and a list of
restriction enzyme names. Each sequence in the input fastq file is in-silico digested with the enzymes specified,
and any resulting sequences beneath a minimum length (default 100 bases) are discarded. All other sequences are
written to an output fastq file (compressed with gzip, by default).
"""

import logging
import sys
import argparse
import os
import socket
import gzip
import re
from Bio import SeqIO, Restriction
import datetime
import itertools
from collections import defaultdict
import textwrap


# f-strings will produce a 'SyntaxError: invalid syntax' error if not supported by Python version:
f'Must be using Python 3.6 or higher.'

if sys.version_info[0:2] < (3, 6):
    sys.exit(f'Must be using Python 3.6 or higher. You are using version {sys.version_info[0]}.{sys.version_info[1]}.')

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")


########################################################################################################################
########################################################################################################################
# Get current working directory and host name:

cwd = os.getcwd()
host = socket.gethostname()


# Configure logger:
def setup_logger(name, log_file, console_level=logging.INFO, file_level=logging.DEBUG,
                 logger_object_level=logging.DEBUG):
    """
    Function to create a logger instance.

    By default, logs level DEBUG and above to file.
    By default, logs level INFO and above to stderr and file.

    :param string name: name for the logger instance
    :param string log_file: filename for log file
    :param string console_level: level for logging to console
    :param string file_level: level for logging to file
    :param string logger_object_level: level for logger object
    :return: a logger object
    """

    # Get date and time string for log filename:
    date_and_time = datetime.datetime.now().strftime("%Y-%m-%d-%H_%M_%S")

    # Log to file:
    file_handler = logging.FileHandler(f'{log_file}_{date_and_time}.log', mode='w')
    file_handler.setLevel(file_level)
    file_format = logging.Formatter('%(asctime)s - %(filename)s - %(name)s - %(funcName)s - %(levelname)s - %('
                                    'message)s')
    file_handler.setFormatter(file_format)

    # Log to Terminal (stdout):
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(console_level)
    console_format = logging.Formatter('%(message)s')
    console_handler.setFormatter(console_format)

    # Setup logger:
    logger_object = logging.getLogger(name)
    logger_object.setLevel(logger_object_level)  # Default level is 'WARNING'

    # Add handlers to the logger
    logger_object.addHandler(console_handler)
    logger_object.addHandler(file_handler)

    return logger_object


# Create logger(s):
logger = setup_logger(__name__, 'digest_then_size_select')


########################################################################################################################
########################################################################################################################
# Define functions:

def createfolder(directory):
    """
    Attempts to create a directory named after the name provided, and provides an error message on failure

    :param string directory: name/path to the directory to be created
    """

    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        logger.info(f'Error: Creating directory: {directory}')


def file_exists_and_not_empty(file_name):
    """
    Check if file exists and is not empty by confirming that its size is not 0 bytes

    :param string file_name:
    """

    # Check if file exist and is not empty
    return os.path.isfile(file_name) and not os.path.getsize(file_name) == 0


def gunzip(file):
    """
    Unzips a .gz file unless unzipped file already exists

    :param string file:
    """

    expected_unzipped_file = re.sub('.gz', '', file)
    if not file_exists_and_not_empty(expected_unzipped_file):
        with open(expected_unzipped_file, 'w') as outfile:
            with gzip.open(file, 'rt') as infile:
                outfile.write(infile.read())
        os.remove(file)

    return expected_unzipped_file


def pairwise(iterable):
    """
    s -> (s0,s1), (s1,s2), (s2, s3), ...
    Used to iterate over overlapping pairs of indices.

    :param list iterable:
    """

    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


def cut_filter_and_report(seq_iterator, verbose_log, enzyme_batch, output_folder_name, filter_fastq_file_name,
                          min_length, uncompressed_output):
    """

    :param seq_iterator:
    :param verbose_log:
    :param enzyme_batch:
    :param output_folder_name:
    :param filter_fastq_file_name:
    :param min_length:
    :param bool uncompressed_output: if True, write an uncompressed output fastq file
    :return:
    """

    # Create a dictionary to count number of sequences cut by each enzyme:
    sequence_enzyme_counter_dict = defaultdict(int)
    for enzyme in enzyme_batch:
        sequence_enzyme_counter_dict[str(enzyme)] = 0

    # Set counters:
    sequence_counter = 0
    cumulative_sequence_counter = 0
    sequence_counter_cut_at_least_once = 0

    if uncompressed_output:
        filtered_fastq_handle = open(f'{output_folder_name}/{filter_fastq_file_name}', 'w')
    else:
        filtered_fastq_handle = gzip.open(f'{output_folder_name}/{filter_fastq_file_name}.gz', 'wt')

    for seq in seq_iterator:
        sequence_counter += 1
        
        if not sequence_counter % 10000:
            sys.stderr.write(f'\r{"[INFO]:":10} Number of sequences processed: {sequence_counter}')

        cut_site_dict = enzyme_batch.search(seq.seq, linear=True)

        for enzyme, cut_site_locations in cut_site_dict.items():
            if cut_site_locations:
                sequence_enzyme_counter_dict[str(enzyme)] += 1

        ordered_sites = sorted({value for values in cut_site_dict.values() for value in values})

        if not ordered_sites:
            SeqIO.write(seq, filtered_fastq_handle, 'fastq')
            # logger.info(f'{"[INFO]:":10} No cut sites found for sequence {seq.name}')
        else:
            sequence_counter_cut_at_least_once += 1
            ordered_sites.insert(0, 1)  # insert a coordinate for the first base of the sequence
            ordered_sites.append(len(seq))  # insert a coordinate for the last base of the sequence
            index_pairs_above_length_cutoff = []

            for index1, index2 in pairwise(ordered_sites):
                if index2 - index1 < min_length:
                    if verbose_log:
                        logger.debug(f'{"[INFO]:":10} Fragment from bases {index1} - {index2} for {seq.name} is '
                                     f'shorter than the minimum length specified ({min_length}). Discarding '
                                     f'fragment...')
                else:
                    index_pairs_above_length_cutoff.append([index1, index2])

            if verbose_log:
                logger.debug(f'{"[INFO]:":10} Retained fragment indexes for {seq.name}:'
                             f' {index_pairs_above_length_cutoff}')

            fragments_to_write = []
            for fragment_index_pair in index_pairs_above_length_cutoff:
                fragment_seq = seq[fragment_index_pair[0] - 1: fragment_index_pair[1] - 1]
                fragments_to_write.append(fragment_seq)

            # Write sequences to file:
            SeqIO.write(fragments_to_write, filtered_fastq_handle, 'fastq')

            if len(fragments_to_write) > 1 and verbose_log:
                logger.debug(f'{"[INFO]:":10} Following in-silico digestion there is more than one fragment '
                             f'above the minimum length threshold for sequence {seq.name}.')

    filtered_fastq_handle.close()
    # logger.info(f'\n')

    # Write stats to screen and file:
    logger.info(f'\n{"[INFO]:":10} Total number of input sequences processed: {sequence_counter}')
    logger.info(f'{"[INFO]:":10} Number of input sequences cut at least once: {sequence_counter_cut_at_least_once}/'
                f'{sequence_counter} ({((sequence_counter_cut_at_least_once / sequence_counter) * 100):.2f}%)')

    for enzyme, sequence_cut_count in sequence_enzyme_counter_dict.items():
        logger.info(f'{"[INFO]:":10} Number of sequences cut with enzyme {enzyme}: {sequence_cut_count}'
                    f'/{sequence_counter} ({((sequence_cut_count / sequence_counter) * 100):.2f}%)')


def digest_and_filter_by_size(fastq_file, enzymes, min_length, verbose_log, uncompressed_output):
    """
    Takes a fastq file, optionally in *.gz format.

    :param str fastq_file: path to a fastq file - can be in *.gz format
    :param list enzymes: list of Bio.Restriction enzyme classes for in-silico digestion
    :param int min_length: minimum number of nucleotides required for a sequence to be retained
    :param bool verbose_log: if True, write very verbose log file
    :param bool uncompressed_output: if True, write an uncompressed output fastq file
    :return:
    """

    # Create output directory:
    output_folder_name = f'digested_and_filtered_fastq_output'
    createfolder(output_folder_name)

    # Get output file name:
    gzipped = False
    file, ext = os.path.splitext(os.path.basename(fastq_file))
    if ext == '.gz':
        gzipped = True
        file, ext = os.path.splitext(file)

    filter_fastq_file_name = f'{file}_filtered{ext}'

    # Create a restriction batch:
    enzyme_batch = Restriction.RestrictionBatch(enzymes)

    if gzipped:
        logger.info(f'{"[INFO]:":10} Processing gzipped fastq input file: {os.path.basename(fastq_file)}')

        with gzip.open(fastq_file, 'rt') as handle:
            seqs = SeqIO.parse(handle, 'fastq')

            cut_filter_and_report(seqs,
                                  verbose_log,
                                  enzyme_batch,
                                  output_folder_name,
                                  filter_fastq_file_name,
                                  min_length,
                                  uncompressed_output)

    else:
        logger.info(f'{"[INFO]:":10} Processing fastq input file {os.path.basename(fastq_file)}')
        seqs = SeqIO.parse(fastq_file, 'fastq')

        cut_filter_and_report(seqs,
                              verbose_log,
                              enzyme_batch,
                              output_folder_name,
                              filter_fastq_file_name,
                              min_length,
                              uncompressed_output)


########################################################################################################################
########################################################################################################################

def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fastq_file',
                        type=str,
                        help='A file of sequences in *.fastq format. Can be gzipped (file name suffix .gz)')
    parser.add_argument('min_length',
                        type=int,
                        default=100,
                        help='Minimum number of bases for a sequence to be retained following in-silico digestion. '
                             'Default is %(default)s')
    parser.add_argument('--restriction_enzymes',
                        nargs='+',
                        required=True,
                        help='One or more restriction enzyme names (e.g. "EcoRI") separated by spaces. Each input DNA '
                             'sequence will be in-silico digested using all the enzymes listed.')
    parser.add_argument('--uncompressed_output',
                        action='store_true',
                        default=False,
                        help='If this flag is used, the output fastq file will be written in uncompressed text format. '
                             'Note that the file size will be MUCH larger then a compressed (suffix .gz) file. A '
                             'compressed file will be written by default.')
    parser.add_argument('--verbose_log',
                        action='store_true',
                        default=False,
                        help='If this flag is used, additional details will be written to the log file for _every_ '
                             'sequence processed. Note that this will produce very large log files and will make '
                             'processing much slower!')
    parser.add_argument('--version', '-v',
                        dest='version',
                        action='version',
                        version='%(prog)s v0.0.1',
                        help='Print the script version number.')

    results = parser.parse_args()
    return results


########################################################################################################################
########################################################################################################################
# Run script:

def main():
    args = parse_arguments()
    logger.info(f'Running {__name__} with: {args}')

    # Check that the input fastq file exists and is not empty:
    if not file_exists_and_not_empty(args.fastq_file):
        sys.exit(f'{"[ERROR]:":10} The fastq file "{args.fastq_file}" either does not exist, or is empty. Please '
                 f'check your input!')

    # Check enzyme(s) can be found, and capture enzyme objects in a list:
    all_enzyme_classes = []

    for enzyme in args.restriction_enzymes:
        try:

            enzyme_class = getattr(Restriction, enzyme)
            all_enzyme_classes.append(enzyme_class)

            fill = textwrap.fill(
                f'{"[INFO]:":10} User-specified enzyme "{enzyme}" found, with cut site {enzyme_class.site}. Blunt:'
                f' {enzyme_class.is_blunt()}, 5-prime overhang: {enzyme_class.is_5overhang()}, 3-prime overhang:'
                f' {enzyme_class.is_3overhang()}. Cut site graphic: {enzyme_class.elucidate()}.',
                width=70, subsequent_indent=' ' * 11)

            logger.info(fill)

        except AttributeError:
            sys.exit(f'{"[ERROR]:":10} No enzyme called "{enzyme}" found! Please check enzyme name, and use correct '
                     f'upper-case Latin numbers (e.g. EcoRI, rather than EcoR1).')

    fill = textwrap.fill(
        f'{"[INFO]:":10} Graphic key: the ^ symbol refers to the position of the cut in the sense strand of the '
        f'sequence, _ to the cut on the anti-sense or complementary strand. ^_ means blunt.',
        width=70, subsequent_indent=' ' * 11)

    logger.info(fill)

    # Digest and filter, and write output files and stats:
    digest_and_filter_by_size(args.fastq_file,
                              all_enzyme_classes,
                              args.min_length,
                              args.verbose_log,
                              args.uncompressed_output)


########################################################################################################################
########################################################################################################################

if __name__ == '__main__':
    if not len(sys.argv) >= 1:
        print(__doc__)
        sys.exit()
    sys.exit(main())

########################################################################################################################
########################################################################################################################
