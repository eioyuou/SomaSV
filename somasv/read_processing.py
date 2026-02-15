import pysam
import numpy as np
from math import floor, ceil

from somasv.cigar_parsing import process_cigar_operations_long, parse_cigar_string
from somasv.breakpoints import build_sv_candidate
from somasv.supplementary import (
    splitreadlist, get_supplementary_breakpoints_new, extract_read_phasing_info
)
from somasv.coverage import (
    add_coverage_to_arrays, create_coverage_arrays, coverage_arrays_to_dict
)


def cigarread(read, potential_breakpoints, start, end, haplotypes, phase_sets, sample, platform='HIFI', min_size=30,
              max_gap=30):
    """
    Process a read's CIGAR string to extract large insertions and deletions.
    """
    chr_name = read.reference_name
    read_name = read.query_name

    del_intervals, ins_intervals = process_cigar_operations_long(read, start, end, platform, min_size, max_gap)

    for del_cigar in del_intervals:
        start_pos, end_pos, length = del_cigar
        potential_breakpoints.setdefault(chr_name, []).append(build_sv_candidate(
            [{'chr': chr_name, 'loc': start_pos}, {'chr': chr_name, 'loc': end_pos}],
            'CIGAR', read_name, read.mapping_quality, sample, "+-", haplotypes, phase_sets
        ))

    for ins_cigar in ins_intervals:
        start_pos, _, length, insert_seq = ins_cigar
        potential_breakpoints.setdefault(chr_name, []).append(build_sv_candidate(
            [{'chr': chr_name, 'loc': start_pos}, {'chr': chr_name, 'loc': start_pos}],
            'CIGAR', read_name, read.mapping_quality, sample, "<INS>", haplotypes, phase_sets, insert_seq
        ))

    return potential_breakpoints


def get_potential_breakpoints(
        aln_filename, length, mapq, sample, contig_order, contig, start, end,
        coverage_binsize, platform, single_bnd=False, single_bnd_min_length=None, single_bnd_max_mapq=None):
    """ Iterate through alignment file, tracking potential breakpoints """

    potential_breakpoints = {}
    aln_file = pysam.AlignmentFile(aln_filename, "rb")
    min_length = max((length - floor(length / 2.5)), 0) if sample == 'normal' else length
    mapq = min((mapq - ceil(mapq / 2)), 0) if sample == 'normal' else mapq
    min_softclip = 20
    max_position = aln_file.get_reference_length(contig) if contig in aln_file.references else 0
    coverage_arrays = create_coverage_arrays(max_position, coverage_binsize)
    read_total_length = 0
    for read in aln_file.fetch(contig, start, end):
        if read.is_unmapped:
            continue
        if read.is_secondary and platform == "HIFI":
            continue

        haplotype, phase_set = extract_read_phasing_info(read)

        add_coverage_to_arrays(coverage_arrays, read.reference_start,
                               read.reference_end, haplotype, coverage_binsize)

        if read.query_length:
            read_total_length += read.query_length
        if read.mapping_quality < mapq and platform == "HIFI":
            continue

        haplotypes = [haplotype, haplotype]
        phase_sets = [phase_set, phase_set]

        cigarread(read, potential_breakpoints, start, end, haplotypes, phase_sets,
                  sample, platform=platform, min_size=min_length)

        if read.is_secondary or read.is_unmapped:
            continue
        if read.mapping_quality < mapq:
            continue
        if not read.is_supplementary and read.has_tag('SA'):
            split_read = splitreadlist(read, contig_order, mapq)
            chimeric_breakpoints = get_supplementary_breakpoints_new(
                read, split_read, contig_order, sample, haplotype, phase_set)

            for bp in chimeric_breakpoints:
                potential_breakpoints.setdefault(bp['start_chr'], []).append(bp)

    aln_file.close()
    del aln_file
    coverage_dict = coverage_arrays_to_dict(coverage_arrays, contig)
    return potential_breakpoints, coverage_dict, read_total_length


def get_potential_breakpoints_task(args):
    try:
        aln_filename, length, mapq, sample, contig_order, contig, start, end, coverage_binsize, platform = args

        print(f"Starting: {sample} - {contig}:{start}-{end}")

        potential_breakpoints, coverage_dict, read_total_length = get_potential_breakpoints(
            aln_filename, length, mapq, sample, contig_order, contig, start, end, coverage_binsize, platform
        )

        print(f"Completed: {sample} - {contig}:{start}-{end}")

        return potential_breakpoints, coverage_dict, read_total_length

    except Exception as e:
        print(f"Error in task {sample} - {contig}:{start}-{end}: {e}")
        return {}


def get_potential_breakpoints_contig(
        aln_filename, length, mapq, sample, contig_order, start, end,
        coverage_binsize, platform, single_bnd=False, single_bnd_min_length=None, single_bnd_max_mapq=None):
    """ Iterate through alignment file for contig analysis """

    potential_breakpoints = {}
    aln_file = pysam.AlignmentFile(aln_filename, "rb")
    min_length = max((length - floor(length / 2.5)), 0) if sample == 'normal' else length
    mapq = min((mapq - ceil(mapq / 2)), 0) if sample == 'normal' else mapq
    min_softclip = 20
    read_total_length = 0
    print(platform)
    for read in aln_file.fetch():
        if read.is_unmapped:
            continue

        haplotype, phase_set = extract_read_phasing_info(read)

        if read.query_length:
            read_total_length += read.query_length

        haplotypes = [0, 0]
        phase_sets = [0, 0]

        cigarread(read, potential_breakpoints, start, end, haplotypes, phase_sets,
                  sample, platform=platform, min_size=min_length)

        if read.is_secondary or read.is_unmapped:
            continue
        if read.mapping_quality < mapq:
            continue
        if not read.is_supplementary and read.has_tag('SA'):
            split_read = splitreadlist(read, contig_order, mapq)
            chimeric_breakpoints = get_supplementary_breakpoints_new(
                read, split_read, contig_order, sample, haplotype, phase_set)

            for bp in chimeric_breakpoints:
                potential_breakpoints.setdefault(bp['start_chr'], []).append(bp)

    aln_file.close()
    del aln_file

    final_potential_breakpoints = []
    for key, value in potential_breakpoints.items():
        if key not in contig_order:
            print(f"Skipping unknown contig: {key}")
            continue
        final_potential_breakpoints.extend(
            sorted(value, key=lambda x: (contig_order.index(x['start_chr']), x['start_loc']))
        )

    return final_potential_breakpoints