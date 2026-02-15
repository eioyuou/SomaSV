import re
from somasv.cigar_parsing import parse_cigar_string
from somasv.breakpoints import build_sv_candidate


def splitreadlist(read, contig_order, mapq):
    sv_list = []
    if not read.has_tag('SA'):
        return []

    rawsalist = read.get_tag('SA').split(';')
    for sa in rawsalist[:-1]:
        sainfo = sa.split(',')
        tmpcontig, tmprefstart, strand, cigar, sup_mapq, nm = sainfo[0], int(sainfo[1]), str(sainfo[2]), sainfo[3], int(
            sainfo[4]), int(sainfo[5])
        if sup_mapq <= mapq:
            continue
        if tmpcontig not in contig_order:
            continue
        is_reverse = True if strand == '-' else False
        read_data = parse_cigar_string(tmpcontig, cigar, tmprefstart, read.query_name, sup_mapq, is_reverse)
        sv_list.append(read_data)

    return sv_list


def create_breakpoint_pair(curr_alignment, next_alignment, softclip_min_mapq, contig_order):
    """
    Create a breakpoint pair between two adjacent alignments
    """
    if curr_alignment['is_reverse']:
        first_bp = {
            'chr': curr_alignment['chr'],
            'loc': curr_alignment['ref_start'],
            'bp_notation': '-',
            'is_reliable_mapping': curr_alignment['mapping_quality'] >= softclip_min_mapq and curr_alignment[
                'chr'] in contig_order,
            'primary': curr_alignment.get('is_primary', False),
            'query_pos_start': curr_alignment['query_start'],
            'query_pos_end': curr_alignment['query_end'],
            'mapping_quality': curr_alignment['mapping_quality']
        }
    else:
        first_bp = {
            'chr': curr_alignment['chr'],
            'loc': curr_alignment['ref_end'],
            'bp_notation': '+',
            'is_reliable_mapping': curr_alignment['mapping_quality'] >= softclip_min_mapq and curr_alignment[
                'chr'] in contig_order,
            'primary': curr_alignment.get('is_primary', False),
            'query_pos_start': curr_alignment['query_start'],
            'query_pos_end': curr_alignment['query_end'],
            'mapping_quality': curr_alignment['mapping_quality']
        }

    if next_alignment['is_reverse']:
        second_bp = {
            'chr': next_alignment['chr'],
            'loc': next_alignment['ref_end'],
            'bp_notation': '+',
            'is_reliable_mapping': next_alignment['mapping_quality'] >= softclip_min_mapq and next_alignment[
                'chr'] in contig_order,
            'primary': next_alignment.get('is_primary', False),
            'query_pos_start': next_alignment['query_start'],
            'query_pos_end': next_alignment['query_end'],
            'mapping_quality': next_alignment['mapping_quality']
        }
    else:
        second_bp = {
            'chr': next_alignment['chr'],
            'loc': next_alignment['ref_start'],
            'bp_notation': '-',
            'is_reliable_mapping': next_alignment['mapping_quality'] >= softclip_min_mapq and next_alignment[
                'chr'] in contig_order,
            'primary': next_alignment.get('is_primary', False),
            'query_pos_start': next_alignment['query_start'],
            'query_pos_end': next_alignment['query_end'],
            'mapping_quality': next_alignment['mapping_quality']
        }

    return [first_bp, second_bp]


def assemble_sv_breakpoints(read, sample, contig_order, breakpoint_pairs, haplotype, phase_set):
    supplementary_breakpoints = []
    index = 0
    buffer = 50

    while index < len(breakpoint_pairs):
        curr_start, curr_end = breakpoint_pairs[index]

        if curr_start['chr'] not in contig_order or curr_end['chr'] not in contig_order:
            index += 1
            continue

        if curr_end['is_reliable_mapping'] and curr_start['is_reliable_mapping']:
            location = [{'chr': curr_start['chr'], 'loc': curr_start['loc']},
                        {'chr': curr_end['chr'], 'loc': curr_end['loc']}]
            mapq = min(curr_start['mapping_quality'], curr_end['mapping_quality'])
            if curr_start['chr'] == curr_end['chr']:
                if curr_start['loc'] <= curr_end['loc']:
                    haplotypes = [haplotype, None] if curr_start['primary'] else [None, haplotype]
                    phase_sets = [phase_set, None] if curr_start['primary'] else [None, phase_set]
                    bp_pattern = curr_start['bp_notation'] + curr_end['bp_notation']

                    if bp_pattern == "-+":
                        ref_span = abs(curr_end['loc'] - curr_start['loc'])
                        if ref_span < buffer:
                            insert_seq = read.query_sequence[curr_start['query_pos_start']:curr_end['query_pos_end']]
                            supplementary_breakpoints.append(build_sv_candidate(
                                [{'chr': curr_start['chr'], 'loc': curr_start['loc']}] * 2,
                                "SUPPLEMENTARY", read.query_name, mapq, sample,
                                '<INS>',
                                haplotypes, phase_sets,
                                insert_seq
                            ))
                        else:
                            supplementary_breakpoints.append(build_sv_candidate(
                                location, "SUPPLEMENTARY", read.query_name, mapq, sample,
                                bp_pattern,
                                haplotypes, phase_sets
                            ))
                    else:
                        supplementary_breakpoints.append(build_sv_candidate(
                            location, "SUPPLEMENTARY", read.query_name, mapq, sample,
                            bp_pattern,
                            haplotypes, phase_sets
                        ))
                else:
                    haplotypes = [haplotype, None] if curr_end['primary'] else [None, haplotype]
                    phase_sets = [phase_set, None] if curr_end['primary'] else [None, phase_set]
                    bp_pattern = curr_end['bp_notation'] + curr_start['bp_notation']

                    if bp_pattern == "-+":
                        ref_span = abs(curr_start['loc'] - curr_end['loc'])
                        if ref_span < buffer:
                            insert_seq = read.query_sequence[curr_end['query_pos_start']:curr_start['query_pos_end']]
                            supplementary_breakpoints.append(build_sv_candidate(
                                [{'chr': curr_end['chr'], 'loc': curr_end['loc']}] * 2,
                                "SUPPLEMENTARY", read.query_name, mapq, sample,
                                '<INS>',
                                haplotypes, phase_sets,
                                insert_seq
                            ))
                        else:
                            supplementary_breakpoints.append(build_sv_candidate(
                                location[::-1], "SUPPLEMENTARY", read.query_name, mapq, sample,
                                bp_pattern,
                                haplotypes, phase_sets
                            ))
                    else:
                        supplementary_breakpoints.append(build_sv_candidate(
                            location[::-1], "SUPPLEMENTARY", read.query_name, mapq, sample,
                            bp_pattern,
                            haplotypes, phase_sets
                        ))
            elif contig_order.index(curr_start['chr']) <= contig_order.index(curr_end['chr']):
                haplotypes = [haplotype, None] if curr_start['primary'] else [None, haplotype]
                phase_sets = [phase_set, None] if curr_start['primary'] else [None, phase_set]
                supplementary_breakpoints.append(build_sv_candidate(
                    location, "SUPPLEMENTARY", read.query_name, mapq, sample,
                    curr_start['bp_notation'] + curr_end['bp_notation'],
                    haplotypes, phase_sets
                ))
            else:
                haplotypes = [haplotype, None] if curr_end['primary'] else [None, haplotype]
                phase_sets = [phase_set, None] if curr_end['primary'] else [None, phase_set]
                supplementary_breakpoints.append(build_sv_candidate(
                    location[::-1], "SUPPLEMENTARY", read.query_name, mapq, sample,
                    curr_end['bp_notation'] + curr_start['bp_notation'],
                    haplotypes, phase_sets
                ))
        index += 1

    return supplementary_breakpoints


def detect_tandem_duplications(all_alignments, read, sample, haplotype, phase_set):
    """Detect tandem duplication events"""
    duplications = []

    for i in range(len(all_alignments)):
        for j in range(i + 1, len(all_alignments)):
            align1 = all_alignments[i]
            align2 = all_alignments[j]

            if (align1['chr'] == align2['chr'] and
                    align1['is_reverse'] == align2['is_reverse']):

                overlap_start = max(align1['ref_start'], align2['ref_start'])
                overlap_end = min(align1['ref_end'], align2['ref_end'])

                if overlap_start < overlap_end:
                    overlap_length = overlap_end - overlap_start
                    if overlap_length >= 30:
                        location = [
                            {'chr': align1['chr'], 'loc': overlap_start},
                            {'chr': align1['chr'], 'loc': overlap_end}
                        ]

                        mapq = min(align1['mapping_quality'], align2['mapping_quality'])

                        haplotypes = [haplotype, None] if align1.get('is_primary', False) else [None, haplotype]
                        phase_sets = [phase_set, None] if align1.get('is_primary', False) else [None, phase_set]

                        bp_pattern = "-+"

                        duplication = build_sv_candidate(
                            location, "SUPPLEMENTARY", read.query_name, mapq, sample,
                            bp_pattern, haplotypes, phase_sets
                        )

                        duplications.append(duplication)

    return duplications


def get_supplementary_breakpoints_new(read, split_read, contig_order, sample, haplotype, phase_set):
    """
    Revised version: correctly handle breakpoints of supplementary alignments
    """
    breakpoint_pairs = []
    potential_breakpoints = []

    chr_name = read.reference_name
    ref_start = read.reference_start
    read_name = read.query_name
    cigarstring = read.cigarstring
    is_reverse = read.is_reverse
    softclip_min_mapq = 5

    primary_read = parse_cigar_string(chr_name, cigarstring, ref_start, read_name, read.mapping_quality, is_reverse)

    split_read = sorted(split_read, key=lambda d: (d['query_start'], d['query_end']))

    all_alignments = split_read + [primary_read]
    all_alignments = sorted(all_alignments, key=lambda d: (d['query_start'], d['query_end']))

    for alignment in all_alignments:
        alignment['is_primary'] = (alignment == primary_read)

    for i in range(len(all_alignments) - 1):
        curr_alignment = all_alignments[i]
        next_alignment = all_alignments[i + 1]

        breakpoint_pair = create_breakpoint_pair(curr_alignment, next_alignment, softclip_min_mapq, contig_order)
        if breakpoint_pair:
            breakpoint_pairs.append(breakpoint_pair)

    duplications = detect_tandem_duplications(all_alignments, read, sample, haplotype, phase_set)
    if duplications:
        potential_breakpoints.extend(duplications)

    potential_breakpoints.extend(
        assemble_sv_breakpoints(read, sample, contig_order, breakpoint_pairs, haplotype, phase_set))

    return potential_breakpoints


def extract_read_phasing_info(read):
    """
    Extract phase (haplotype) and phase set information.
    """
    haplotype, phase_set = None, None
    if read.is_supplementary or read.is_secondary:
        return haplotype, phase_set

    if read.has_tag('HP'):
        haplotype = read.get_tag('HP')

    if read.has_tag('PS'):
        phase_set = read.get_tag('PS')

    return haplotype, phase_set