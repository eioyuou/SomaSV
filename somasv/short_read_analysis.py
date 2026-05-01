import re
import pysam
import numpy as np
from somasv.cigar_parsing import parse_cigar_string




SVET_UNKNOWN = 0
SVET_CIGAR = 1
SVET_PAIR = 2
SVET_SPLIT_ALIGN = 3
SVET_SEMIALIGN = 4
SVET_SHADOW = 5
SVET_LOCAL_PAIR = 6
SVET_SIZE = 7


def create_sv_breakend_short(parsed_cigar_result, state='SVBEND_UNKNOWN'):
    breakend = {}
    breakend["interval"] = {
        "chr": parsed_cigar_result["chr"],
        "range": {
            "begin_pos": parsed_cigar_result["ref_start"],
            "end_pos": parsed_cigar_result["ref_end"]
        }
    }
    breakend["state"] = state
    breakend["lowres_evidence"] = [0] * SVET_SIZE
    return breakend


def cigarstring_to_path_short(cigarstring):
    return [{"type": op, "length": int(length_str)}
            for length_str, op in re.findall(r'(\d+)([MIDNSHP=X])', cigarstring)]


def is_segment_align_match_short(segment_type):
    return segment_type in ['M', '=', 'X']


def is_segment_type_indel_short(segment_type):
    return segment_type in ['I', 'D']


def is_segment_type_read_length_short(segment_type):
    return segment_type in ['M', 'I', 'S', '=', 'X']


def is_segment_type_ref_length_short(segment_type):
    return segment_type in ['M', 'D', 'N', '=', 'X']


def is_base_match_for_poor_alignment_test_short(a, b):
    if a == 'N' or b == 'N':
        return True
    return a.upper() == b.upper()


def get_scale_short(scaler, size):
    """获取缩放因子"""
    if size <= scaler["min_size"]:
        return 0.0
    elif size >= scaler["max_size"]:
        return 1.0
    else:
        return (size - scaler["min_size"]) / (scaler["max_size"] - scaler["min_size"])


def parse_cigar_string_short(chr_name, cigarstring, ref_start, read_name, mapping_quality, is_reverse=False):
    """
    Parse CIGAR string, calculate reference genome alignment start, end position and query sequence position.
    """
    # Initialize various length information
    result = {
        'chr': chr_name,
        'softclip_left': 0,
        'softclip_right': 0,
        'match_length': 0,
        'insert_length': 0,
        'delete_length': 0,
        'mismatch_length': 0,
        'ref_start': ref_start,
        'ref_end': ref_start,
        'query_start': 0,
        'query_end': 0,
        'read_name': read_name,
        'is_reverse': is_reverse,
        'mapping_quality': mapping_quality
    }

    # Regular expression to match CIGAR operations (e.g.: 5S, 10M, 3I, 2D)
    cigar_ops = re.findall(r'(\d+)([MIDNSHP=X])', cigarstring)

    # For reverse aligned reads, CIGAR operations need to be reversed
    # if is_reverse:
    #     cigar_ops = cigar_ops[::-1]

    query_pos = 0  # Position on query sequence (0-based)
    ref_pos = ref_start  # Current position on reference genome

    # Iterate through CIGAR operations, process each operation's length and type
    for i, (length_str, op) in enumerate(cigar_ops):
        length = int(length_str)  # Convert operation length to integer

        if op == 'M':  # M (match or mismatch)
            result['match_length'] += length
            query_pos += length
            ref_pos += length
        elif op == 'I':  # I (insertion)
            result['insert_length'] += length
            query_pos += length
        elif op == 'D':  # D (deletion)
            result['delete_length'] += length
            ref_pos += length
        elif op == 'S':  # S (softclip)
            if i == 0:  # First S operation is left-side softclip
                result['softclip_left'] = length
            else:  # Last S operation is right-side softclip
                result['softclip_right'] = length
            # Softclip only affects query sequence position, not reference sequence
            query_pos += length
        elif op == '=':  # = (exact match)
            result['match_length'] += length
            query_pos += length
            ref_pos += length
        elif op == 'X':  # X (mismatch)
            result['mismatch_length'] += length
            query_pos += length
            ref_pos += length

        elif op == 'H':  # Hard clip
            # Hard clip不影响坐标，跳过
            pass

        elif op == 'N':  # Skipped region (intron)
            ref_pos += length
            # query_pos不变

    # Calculate ref_end, final reference genome end position
    result['ref_end'] = ref_pos

    # Calculate query sequence start and end positions (1-based)
    result['query_start'] = result['softclip_left']  # Start position is the first position after softclip (1-based)
    result['query_end'] = query_pos - result['softclip_right']  # Final position after removing right softclip

    return result


def is_match_tuple_op_short(op):
    return op in (0, 7, 8)  # M, =, X

def get_flank_lengths_short(cigartuples, inx):
    left = 0
    for j in range(inx - 1, -1, -1):
        op, l = cigartuples[j]
        if is_match_tuple_op_short(op):
            left += l
        else:
            break
    right = 0
    for j in range(inx + 1, len(cigartuples)):
        op, l = cigartuples[j]
        if is_match_tuple_op_short(op):
            right += l
        else:
            break
    return left, right


def process_cigar_operations_short(read, start, end, min_size=40, max_gap=150):
    """
    Process the CIGAR string from a read, identifying large insertions and deletions.
    Returns two lists: deletions and insertions, where insertions include the actual sequence.
    """
    del_intervals = []
    ins_intervals = []
    reference_pos = read.reference_start
    query_pos = 0  # Position in query sequence
    left_flank = 0  # Temporary variable to store length of current operation
    right_flank = 0  # Temporary variable to store length of current operation
    flanking_thread = 5

    for inx, (operation, length) in enumerate(read.cigartuples):
        if operation in [0, 7, 8]:  # Match or alignment match
            reference_pos += length
            query_pos += length
        elif operation == 2:  # Deletion
            if length >= min_size and start <= reference_pos + length and end >= reference_pos:
                left_flank, right_flank = get_flank_lengths_short(read.cigartuples, inx)
                if left_flank < flanking_thread or right_flank < flanking_thread:
                    continue
                del_intervals.append([reference_pos, reference_pos + length, length])
            reference_pos += length
        elif operation == 1:  # Insertion
            if length >= min_size and start <= reference_pos and end >= reference_pos:
                # Extract inserted sequence (from query sequence)
                left_flank, right_flank = get_flank_lengths_short(read.cigartuples, inx)
                if left_flank < flanking_thread or right_flank < flanking_thread:
                    continue
                insert_seq = read.query_sequence[query_pos:query_pos + length]
                ins_intervals.append([reference_pos, reference_pos, length, insert_seq])
            query_pos += length  # Insertion only advances query sequence, reference stays the same

    # Merge intervals
    del_intervals = merge_intervals_short(del_intervals, max_gap=max_gap)
    ins_intervals = merge_intervals_short(ins_intervals, max_gap=max_gap, is_insertion=True)

    return del_intervals, ins_intervals


def merge_intervals_short(intervals, max_gap=150, is_insertion=False):
    """
    Merge intervals that are close to each other within `max_gap`.
    If the intervals are insertions, their sequences will be merged as well.
    """
    merged = []
    if not intervals:
        return merged
    # Sort intervals by start position
    intervals.sort(key=lambda x: x[0])
    # Initialize the first interval
    if is_insertion:
        current_start, current_end, current_length, current_sequence = intervals[0]
    else:
        current_start, current_end, current_length = intervals[0]
    for interval in intervals[1:]:
        if is_insertion:
            start, end, length, sequence = interval
        else:
            start, end, length = interval

        # If intervals overlap or are within max_gap, merge them
        if start - current_end <= max_gap:
            current_end = max(current_end, end)
            current_length += length
            if is_insertion:
                current_sequence += sequence  # Merge inserted sequences
        else:
            # Save the previous interval and start a new one
            if is_insertion:
                merged.append([current_start, current_end, current_length, current_sequence])
                current_start, current_end, current_length, current_sequence = start, end, length, sequence
            else:
                merged.append([current_start, current_end, current_length])
                current_start, current_end, current_length = start, end, length

    # Append the last interval
    if is_insertion:
        merged.append([current_start, current_end, current_length, current_sequence])
    else:
        merged.append([current_start, current_end, current_length])

    return merged


def get_split_sv_candidate_short(primary_read, ref_chr_length, sample_name, start_pos, end_pos,
                           before_breakend_offset, after_breakend_offset,
                           sv_evidence_source, sv_evidence_type, is_complex=False):
    sv = {
        "bp1": create_sv_breakend_short(primary_read),
        "bp2": create_sv_breakend_short(primary_read),
        "sv_evidence_type": sv_evidence_type,
        'sample_name': sample_name,
    }
    local_breakend = sv["bp1"]
    remote_breakend = sv["bp2"]
    ref_length = ref_chr_length[primary_read['chr']]
    local_breakend["interval"]["chr"] = primary_read['chr']
    remote_breakend["interval"]["chr"] = primary_read['chr']
    local_breakend["lowres_evidence"][sv_evidence_source] += 1

    sv["sv_evidence_type"] = sv_evidence_type
    if not is_complex:
        remote_breakend["lowres_evidence"][sv_evidence_source] += 1
        local_breakend["state"] = 'SVBEND_RIGHT_OPEN'
        remote_breakend["state"] = 'SVBEND_LEFT_OPEN'
    else:
        local_breakend["state"] = 'SVBEND_COMPLEX'
        remote_breakend["state"] = 'SVBEND_UNKNOWN'
    beforeBreakend = before_breakend_offset
    afterBreakend = after_breakend_offset
    local_breakend["interval"]["range"]["begin_pos"] = max(0, start_pos - beforeBreakend)
    if not is_complex:
        local_breakend["interval"]["range"]["end_pos"] = min(ref_length, start_pos + afterBreakend)
    else:
        local_breakend["interval"]["range"]["end_pos"] = min(ref_length, end_pos + afterBreakend)
    remote_breakend["interval"]["range"]["begin_pos"] = max(0, end_pos - beforeBreakend)
    remote_breakend["interval"]["range"]["end_pos"] = min(ref_length, end_pos + afterBreakend)

    return sv


def get_sv_candidates_from_read_indels_short(read, ref_chr_length, sample_name, start, end, min_size=50, max_gap=150,
                                       before_breakend_offset=100, after_breakend_offset=100):
    candidate = []
    sv_source = SVET_CIGAR
    chr_name = read.reference_name
    ref_start = read.reference_start
    read_name = read.query_name
    cigarstring = read.cigarstring
    is_reverse = read.is_reverse
    is_read2 = read.is_paired and not read.is_read1
    frag_source = 'read2' if is_read2 else 'read1'
    # Parse primary alignment CIGAR string

    primary_read = parse_cigar_string_short(chr_name, cigarstring, ref_start, read_name, read.mapping_quality, is_reverse)
    del_intervals, ins_intervals = process_cigar_operations_short(read, start, end, min_size=min_size, max_gap=max_gap)

    # Add deletion events to potential_breakpoints
    for del_cigar in del_intervals:
        start_pos, end_pos, length = del_cigar
        sv_evidence_type = 'DEL'
        sv_candidates = get_split_sv_candidate_short(primary_read, ref_chr_length, sample_name, start_pos, end_pos,
                                               before_breakend_offset, after_breakend_offset, sv_source,
                                               sv_evidence_type,
                                               is_complex=False)
        candidate.append(sv_candidates)
    # Add insertion events to potential_breakpoints (including actual inserted sequence)
    for ins_cigar in ins_intervals:
        start_pos, _, length, insert_seq = ins_cigar
        sv_evidence_type = 'INS'
        sv_candidates = get_split_sv_candidate_short(primary_read, ref_chr_length, sample_name, start_pos, start_pos,
                                               before_breakend_offset, after_breakend_offset,
                                               sv_source, sv_evidence_type, is_complex=False)
        sv_candidates['insert_seq'] = insert_seq
        sv_candidates['sample_name'] = sample_name
        candidate.append(sv_candidates)

    return candidate


def is_split_open_downstream_short(readsig):

    lead, trail = readsig['softclip_left'], readsig['softclip_right']
    return lead < trail


def validate_sa_alignment_short(primary_read, sa_alignment):
    """
    Validate SA (Supplementary Alignment) consistency, inspired by Manta-style logic.

    Key criteria:
    - The two alignments (primary & SA) should be on the same chromosome with a reasonable reference distance.
    - Their coverage on the read (query sequence) should generally not exceed expected length.
    - At least one end should display a significant soft-clip, indicating a potential split.
    """

    # Parse primary alignment details
    primary_align = parse_cigar_string_short(
        primary_read.reference_name,
        primary_read.cigarstring,
        primary_read.reference_start,
        primary_read.query_name,
        primary_read.mapping_quality,
        primary_read.is_reverse
    )

    # 1. Chromosome and reference distance check
    if primary_align["chr"] == sa_alignment["chr"]:
        distance = abs(primary_align["ref_start"] - sa_alignment["ref_start"])
        # If too close, it might be noise; too far may indicate incorrect pairing
        if distance < 20:
            return False

    # 2. Query span consistency check
    primary_query_len = primary_align["query_end"] - primary_align["query_start"]
    sa_query_len = sa_alignment["query_end"] - sa_alignment["query_start"]
    total_query_len = primary_query_len + sa_query_len

    # Allow slight overlap, but total query-covered length shouldn't exceed read length significantly
    if total_query_len > primary_read.query_length + 10:
        return False

    # 3. Soft-clipping requirement
    primary_left_clip = primary_align["softclip_left"]
    primary_right_clip = primary_align["softclip_right"]
    sa_left_clip = sa_alignment["softclip_left"]
    sa_right_clip = sa_alignment["softclip_right"]

    # A valid split-read alignment should show significant soft-clip on at least one side
    if (primary_left_clip < 10 and primary_right_clip < 10 and
            sa_left_clip < 10 and sa_right_clip < 10):
        return False

    return True
def validate_sa_alignment_short(primary_read, sa_alignment,
                          min_aligned_bases=20,       # Minimum number of aligned (M/=X) bases per segment
                          max_query_overlap=20,       # Max allowed query overlap between the two parts
                          max_query_gap=20,           # Max allowed query gap between the two parts
                          junction_slop=10,           # Allowed margin around junction in query space
                          max_total_len_excess=30,    # Allowed excess in total query-covered length over read
                          min_softclip_any=8,         # Require soft clipping on either side for nearby events
                          far_ref_distance=1000       # Relax soft-clip check for far or interchromosomal joins
                          ):
    """
    Robust validation of supplementary alignment (SA tag):
    - Based on query contiguity (real split-read behavior), not just reference distance.
    - Both segments must have sufficient aligned bases (matched/mismatched).
    - For nearby rearrangements, soft-clipping at junction is usually expected; relaxed for distant/interchromosomal.
    - Avoids false split-reads where query-covered length significantly exceeds the original read.
    """

    # Parse primary alignment using supporting utility (assumed implemented)
    primary_align = parse_cigar_string_short(
        primary_read.reference_name,
        primary_read.cigarstring,
        primary_read.reference_start,
        primary_read.query_name,
        primary_read.mapping_quality,
        primary_read.is_reverse
    )

    # 1. Check that each part has enough aligned bases (M/=/X count)
    prim_aln_bases = primary_align.get("match_length", 0) + primary_align.get("mismatch_length", 0)
    sa_aln_bases   = sa_alignment.get("match_length", 0) + sa_alignment.get("mismatch_length", 0)
    if prim_aln_bases < min_aligned_bases or sa_aln_bases < min_aligned_bases:
        return False

    # 2. Query adjacency / splicing consistency
    p_qs, p_qe = primary_align["query_start"], primary_align["query_end"]
    s_qs, s_qe = sa_alignment["query_start"], sa_alignment["query_end"]

    # Two junction modes: primary right --> SA left OR SA right --> primary left
    j1 = abs(p_qe - s_qs)
    j2 = abs(s_qe - p_qs)
    junction_ok = (j1 <= junction_slop) or (j2 <= junction_slop)

    # If not tightly joined, fall back on overlap/gap limits
    if not junction_ok:
        overlap = min(p_qe, s_qe) - max(p_qs, s_qs)  # >0 => overlap
        gap = max(p_qs, s_qs) - min(p_qe, s_qe)      # >0 => gap
        if overlap > 0 and overlap > max_query_overlap:
            return False
        if gap > 0 and gap > max_query_gap:
            return False

    # 3. Coverage length consistency: total mapped query length shouldn't drastically exceed read length
    read_len = primary_read.query_length or (max(p_qe, s_qe) - min(p_qs, s_qs))  # fallback if missing
    total_covered = (p_qe - p_qs) + (s_qe - s_qs)
    if total_covered > read_len + max_total_len_excess:
        return False

    # 4. Soft-clipping validation for nearby events
    same_chr = (primary_align["chr"] == sa_alignment["chr"])
    ref_dist = 0
    if same_chr:
        # Use closest reference edge-to-edge distance
        p_rs, p_re = primary_align["ref_start"], primary_align["ref_end"]
        s_rs, s_re = sa_alignment["ref_start"], sa_alignment["ref_end"]
        ref_dist = max(0, max(p_rs, s_rs) - min(p_re, s_re))  # 0 if overlapping

    need_softclip = (not same_chr) and False  # Don't enforce softclip for inter-chromosomal events
    if same_chr and ref_dist <= far_ref_distance:
        need_softclip = True  # For nearby events, typically we expect soft clipping at junction

    if need_softclip:
        # Evaluate both junction orientation patterns:
        p_sc_left, p_sc_right = primary_align["softclip_left"], primary_align["softclip_right"]
        s_sc_left, s_sc_right = sa_alignment["softclip_left"],  sa_alignment["softclip_right"]

        sc_ok = False
        if j1 <= j2:
            # primary right joins SA left
            if max(p_sc_right, s_sc_left) >= min_softclip_any:
                sc_ok = True
        else:
            # SA right joins primary left
            if max(s_sc_right, p_sc_left) >= min_softclip_any:
                sc_ok = True

        # Fallback - allow any side to have clear soft-clipping
        if not sc_ok and max(p_sc_left, p_sc_right, s_sc_left, s_sc_right) >= min_softclip_any:
            sc_ok = True

        if not sc_ok:
            return False

    return True

def splitreadlist_short(read, chrom_list, min_mapq):
    sv_list = []
    if not read.has_tag('SA'):
        return []
    rawsalist = read.get_tag('SA').split(';')[:-1]
    # if len(rawsalist) >=5 :
    #     return []
    for sa in rawsalist[:]:
        sainfo = sa.split(',')
        tmpcontig, tmprefstart, strand, cigar, sup_mapq, nm = sainfo[0], int(sainfo[1]), str(sainfo[2]), sainfo[3], int(
            sainfo[4]), int(sainfo[5])
        if tmpcontig not in chrom_list:
            continue
        tmprefstart -= 1  # Convert to 0-based index
        is_reverse = True if strand == '-' else False
        if sup_mapq < min_mapq:
            continue
        read_data = parse_cigar_string_short(tmpcontig, cigar, tmprefstart, read.query_name, read.mapping_quality, is_reverse)
        if validate_sa_alignment_short(read, read_data):
            sv_list.append(read_data)
    return sv_list


def update_sa_breakend_short(breakend, local_align, ref_chr_length, before_breakend_offset, after_breakend_offset):
    is_split_downstream = is_split_open_downstream_short(local_align)

    if is_split_downstream:
        breakend["state"] = 'SVBEND_RIGHT_OPEN'
        pos = local_align['ref_end']
    else:
        breakend["state"] = 'SVBEND_LEFT_OPEN'
        pos = local_align['ref_start']

    breakend["interval"]["chr"] = local_align["chr"]
    chrom_length = ref_chr_length[local_align["chr"]]

    begin_pos = max(0, pos - before_breakend_offset)
    end_pos = min(chrom_length, pos + after_breakend_offset)

    breakend["interval"]["range"]["begin_pos"] = begin_pos
    breakend["interval"]["range"]["end_pos"] = end_pos


def get_split_sa_candidate_short(local_align, sa_align,
                           ref_chr_length, sample_name, before_breakend_offset, after_breakend_offset):

    sv = {
        "bp1": create_sv_breakend_short(local_align),
        "bp2": create_sv_breakend_short(sa_align),
        "sv_evidence_type": 'SPLIT',
        "sample_name": sample_name,
    }

    local_breakend = sv["bp1"]
    remote_breakend = sv["bp2"]


    local_breakend["lowres_evidence"][SVET_SPLIT_ALIGN] += 1


    update_sa_breakend_short(local_breakend, local_align, ref_chr_length,
                       before_breakend_offset, after_breakend_offset)
    update_sa_breakend_short(remote_breakend, sa_align, ref_chr_length,
                       before_breakend_offset, after_breakend_offset)

    return sv


def extract_split_reads_short(read, ref_chr_length, sample_name, chrom_list, min_mapq=20, before_breakend_offset=100,
                        after_breakend_offset=100,
                        is_transcript_strand_known=False):

    candidates = []

    if read.has_tag('SA'):


        is_read2 = read.is_paired and not read.is_read1
        frag_source = 'read2' if is_read2 else 'read1'
        # Primary alignment CIGAR string information

        chr_name = read.reference_name
        ref_start = read.reference_start
        read_name = read.query_name
        cigarstring = read.cigarstring
        is_reverse = read.is_reverse
        # Parse primary alignment CIGAR string
        primary_read = parse_cigar_string_short(chr_name, cigarstring, ref_start, read_name, read.mapping_quality, is_reverse)

        splitreadlis = splitreadlist_short(read, chrom_list, min_mapq)
        #print(f'len(splitreadlis): {len(splitreadlis)}')

        sv = 0
        for sample in splitreadlis:
            sv = get_split_sa_candidate_short(primary_read, sample, ref_chr_length, sample_name,
                                        before_breakend_offset, after_breakend_offset)

            #sv = get_split_sa_candidate_short(dopt, read, local_align, ral, frag_source, bam_header)
            candidates.append(sv)
    # print(f'candidates: {candidates}')
    return candidates


def is_overlapping_pair_short(read):
    if (not read.is_paired or
            read.is_unmapped or
            read.mate_is_unmapped or
            read.reference_id != read.next_reference_id):
        return False
    read_end = read.reference_end
    mate_start = read.next_reference_start
    mate_end = None
    if read.has_tag('MC'):
        mate_cigar = read.get_tag('MC')
        mate_cigar = re.findall(r'(\d+)([MIDNSHP=X])', mate_cigar)
        mate_ref_len = sum(int(length) for length, op in mate_cigar
                           if op in ['M', 'D', 'N', '=', 'X'])
        mate_end = mate_start + mate_ref_len

    if mate_end is not None:
        return (read.reference_start < mate_end and mate_start < read_end)
    else:
        return (read.reference_start <= mate_start and mate_start < read_end) or \
            (mate_start <= read.reference_start and read.reference_start < mate_start + read.query_length)


def is_read_filtered_short(read, min_mapq=20):

    return (read.is_unmapped or
            read.is_duplicate or
            read.is_qcfail or
            read.is_secondary or
            read.is_supplementary or
            read.mapping_quality < min_mapq)


def is_discordant_pair_short(read, min_isize=50, max_isize=100000):
    """
    Determine whether a read pair is discordant (e.g., suggesting structural abnormality).

    Args:
        read: pysam.AlignedSegment object
        min_isize: minimum expected insert size
        max_isize: maximum expected insert size (e.g., for somatic SVs, can be large, like 100kb)
    """
    # Basic conditions
    if (not read.is_paired or
            read.is_unmapped or
            read.mate_is_unmapped or
            read.is_secondary or
            read.is_supplementary):
        return False

    # Inter-chromosomal pair — clearly discordant
    if read.reference_id != read.next_reference_id:
        return True

    # Same chromosome: check insert size
    insert_size = abs(read.template_length)

    # Insert size too small or too big
    if insert_size < min_isize or insert_size > max_isize:
        return True

    # Orientation check
    # Expected:
    # - FR (forward-reverse): read1 forward, read2 reverse, read1 before read2
    # - RF (reverse-forward): read1 reverse, read2 forward, read1 after read2

    # Discordant if both mates have same orientation (both fwd or both rev)
    if read.is_reverse == read.mate_is_reverse:
        return True

    read_start = read.reference_start
    mate_start = read.next_reference_start

    if read.is_read1:
        # If this is read1
        if not read.is_reverse and read.mate_is_reverse:
            # Expected FR: read1 is forward, mate is reverse, read1 should be before mate
            if read_start >= mate_start:
                return True
        elif read.is_reverse and not read.mate_is_reverse:
            # Expected RF: read1 is reverse, mate is forward, read1 should be after mate
            if read_start <= mate_start:
                return True
    else:
        # If this is read2 (logic is flipped)
        if not read.is_reverse and read.mate_is_reverse:
            # Read2 forward, mate (read1) reverse — should be after mate
            if read_start <= mate_start:
                return True
        elif read.is_reverse and not read.mate_is_reverse:
            # Read2 reverse, mate (read1) forward — should be before mate
            if read_start >= mate_start:
                return True

    return False


def get_pair_orientation_short(read):

    if read.is_reverse and read.mate_is_reverse:
        return "RR"
    elif not read.is_reverse and not read.mate_is_reverse:
        return "FF"
    elif not read.is_reverse and read.mate_is_reverse:
        if read.reference_start < read.next_reference_start:
            return "FR"
        else:
            return "RF"
    else:  # read.is_reverse and not read.mate_is_reverse
        if read.reference_start > read.next_reference_start:
            return "RF"
        else:
            return "FR"


def infer_sv_type_from_discordant_pair_short(read):


    if read.reference_id != read.next_reference_id:
        return "TRA"


    insert_size = abs(read.template_length)
    orientation = get_pair_orientation_short(read)


    if insert_size > 1000:
        return "DEL"

    # 方向异常
    if orientation in ["RR", "FF"]:
        return "INV"
    elif orientation == "RF" and read.reference_start < read.next_reference_start:
        return "DUP"

    return "UNK"




def get_discordant_pair_candidate_short(read, ref_chr_length, sample_name, ref_chr_id_to_name,
                                  before_breakend_offset=500, after_breakend_offset=500):

    # 解析read信息
    chr_name = read.reference_name
    ref_start = read.reference_start
    read_name = read.query_name
    cigarstring = read.cigarstring
    is_reverse = read.is_reverse

    primary_read = parse_cigar_string_short(chr_name, cigarstring, ref_start, read_name,
                                      read.mapping_quality, is_reverse)

    sv_type = infer_sv_type_from_discordant_pair_short(read)

    sv = {
        "bp1": create_sv_breakend_short(primary_read),
        "bp2": create_sv_breakend_short(primary_read),
        "sv_evidence_type": sv_type,
        "sample_name": sample_name,
    }

    local_breakend = sv["bp1"]
    remote_breakend = sv["bp2"]


    local_breakend["interval"]["chr"] = read.reference_name
    local_breakend["lowres_evidence"][SVET_PAIR] += 1


    if read.is_reverse:
        local_breakend["state"] = 'SVBEND_LEFT_OPEN'
        local_pos = read.reference_start
    else:
        local_breakend["state"] = 'SVBEND_RIGHT_OPEN'
        local_pos = read.reference_end

    mate_chr = ref_chr_id_to_name[read.next_reference_id] if read.next_reference_id >= 0 else read.reference_name
    remote_breakend["interval"]["chr"] = mate_chr


    if read.mate_is_reverse:
        remote_breakend["state"] = 'SVBEND_LEFT_OPEN'
        remote_pos = read.next_reference_start
    else:
        remote_breakend["state"] = 'SVBEND_RIGHT_OPEN'
        mate_end = read.next_reference_start + read.query_length  # 粗略估算
        if read.has_tag('MC'):
            mate_cigar = read.get_tag('MC')
            mate_cigar_ops = re.findall(r'(\d+)([MIDNSHP=X])', mate_cigar)
            mate_ref_len = sum(int(length) for length, op in mate_cigar_ops
                               if op in ['M', 'D', 'N', '=', 'X'])
            mate_end = read.next_reference_start + mate_ref_len
        remote_pos = mate_end


    local_ref_length = ref_chr_length[read.reference_name]
    local_breakend["interval"]["range"]["begin_pos"] = max(0, local_pos - before_breakend_offset)
    local_breakend["interval"]["range"]["end_pos"] = min(local_ref_length, local_pos + after_breakend_offset)

    if mate_chr in ref_chr_length:
        remote_ref_length = ref_chr_length[mate_chr]
        remote_breakend["interval"]["range"]["begin_pos"] = max(0, remote_pos - before_breakend_offset)
        remote_breakend["interval"]["range"]["end_pos"] = min(remote_ref_length, remote_pos + after_breakend_offset)
    else:

        remote_breakend["interval"]["range"]["begin_pos"] = max(0, remote_pos - before_breakend_offset)
        remote_breakend["interval"]["range"]["end_pos"] = remote_pos + after_breakend_offset


    sv["insert_size"] = abs(read.template_length) if read.template_length else 0
    sv["orientation"] = get_pair_orientation_short(read)

    return sv


def extract_discordant_pairs_short(read, ref_chr_length, sample_name, ref_chr_id_to_name,
                             min_isize=500, max_isize=10000,
                             before_breakend_offset=500, after_breakend_offset=500):

    candidates = []

    if not (read.is_paired and read.is_read1 and not read.mate_is_unmapped):
        return candidates

    if not is_discordant_pair_short(read, min_isize, max_isize):
        return candidates



    candidate = get_discordant_pair_candidate_short(read, ref_chr_length, sample_name, ref_chr_id_to_name,
                                              before_breakend_offset, after_breakend_offset)
    candidates.append(candidate)

    return candidates





def leading_edge_poor_alignment_length_short_short(align, read_seq, ref_file, contiguous_match_count):

    assert contiguous_match_count > 0

    read_index = 0
    ref_index = align["ref_start"]
    match_length = 0

    try:
        ref = pysam.FastaFile(ref_file)

        ref_start = max(0, align["ref_start"])
        ref_end = min(ref.get_reference_length(align["chr"]), align["ref_end"] + 100)
        ref_seq = ref.fetch(align["chr"], ref_start, ref_end).upper()
        ref_offset = ref_start
    except Exception as e:
        logger.debug(f"Error fetching reference: {e}")
        return 0, align["ref_start"]

    for segment in align["cigarstring"]:
        seg_type = segment["type"]
        seg_length = segment["length"]

        if is_segment_align_match_short(seg_type):
            for seg_pos in range(seg_length):
                read_pos = read_index + seg_pos
                ref_pos = ref_index + seg_pos


                if (read_pos >= len(read_seq) or
                        ref_pos < ref_offset or
                        ref_pos - ref_offset >= len(ref_seq)):
                    break

                read_base = read_seq[read_pos].upper()
                ref_base = ref_seq[ref_pos - ref_offset]

                if is_base_match_for_poor_alignment_test_short(read_base, ref_base):
                    match_length += 1
                    if match_length >= contiguous_match_count:
                        leading_length = read_pos - (match_length - 1)
                        leading_ref_pos = ref_pos - (match_length - 1)
                        return leading_length, leading_ref_pos
                else:
                    match_length = 0

        elif is_segment_type_indel_short(seg_type):
            match_length = 0


        if is_segment_type_read_length_short(seg_type):
            read_index += seg_length
        if is_segment_type_ref_length_short(seg_type):
            ref_index += seg_length

    return 0, align["ref_start"]


def trailing_edge_poor_alignment_length_short_short(align, read_seq, ref_file, contiguous_match_count):
    assert contiguous_match_count != 0
    read_size = len(read_seq)
    read_index = read_size - 1
    ref = pysam.FastaFile(ref_file)
    ref_seq = ref.fetch(align["chr"], align["ref_start"], align["ref_end"] + 20)
    ref_pos_offset = align["ref_start"]
    path_ref_length = align["ref_end"] - align["ref_start"]
    ref_index = ref_pos_offset + path_ref_length - 1
    match_length = 0
    for segment in reversed(align["cigarstring"]):
        seg_type = segment["type"]
        seg_length = segment["length"]
        if is_segment_align_match_short(seg_type):
            for seg_pos in range(seg_length):
                if (read_index - seg_pos >= 0 and
                        ref_index - seg_pos >= ref_pos_offset and
                        ref_index - seg_pos - ref_pos_offset < len(ref_seq)):
                    read_pos = read_index - seg_pos
                    ref_pos = ref_index - seg_pos - ref_pos_offset
                    if read_pos < len(read_seq) and ref_pos < len(ref_seq):
                        if is_base_match_for_poor_alignment_test_short(read_seq[read_pos], ref_seq[ref_pos]):
                            match_length += 1
                            if match_length >= contiguous_match_count:
                                trailing_length = (read_size - (read_index - seg_pos)) - match_length
                                trailing_ref_pos = (ref_index - seg_pos) + match_length
                                return trailing_length, trailing_ref_pos
                        else:
                            match_length = 0
        elif is_segment_type_indel_short(seg_type):
            match_length = 0
        if is_segment_type_read_length_short(seg_type):
            read_index -= seg_length
        if is_segment_type_ref_length_short(seg_type):
            ref_index -= seg_length
    trailing_length = 0
    trailing_ref_pos = align["ref_end"]
    return trailing_length, trailing_ref_pos


def matchify_edge_soft_clip_short(align):

    matched_align = dict(align)
    matched_align["cigarstring"] = align["cigarstring"].copy()


    if (matched_align["cigarstring"] and
            matched_align["cigarstring"][0]["type"] == 'S'):
        clip_length = matched_align["cigarstring"][0]["length"]

        matched_align["ref_start"] -= clip_length

        matched_align["cigarstring"][0]["type"] = 'M'


    if (matched_align["cigarstring"] and
            matched_align["cigarstring"][-1]["type"] == 'S'):
        # 将S转换为M，ref_end会在重新计算时自动扩展
        matched_align["cigarstring"][-1]["type"] = 'M'


    ref_pos = matched_align["ref_start"]
    for segment in matched_align["cigarstring"]:
        if segment["type"] in ['M', 'D', 'N', '=', 'X']:
            ref_pos += segment["length"]
    matched_align["ref_end"] = ref_pos

    return matched_align


def edge_poor_alignment_length_short(align, read_seq, ref_file, contiguous_match_count):

    leading_length, leading_ref_pos = 0, align["ref_start"]
    trailing_length, trailing_ref_pos = 0, align["ref_end"]
    if align['softclip_left'] > 0:
        leading_length, leading_ref_pos = leading_edge_poor_alignment_length_short_short(
            align, read_seq, ref_file, contiguous_match_count)
    if align['softclip_right'] > 0:
        trailing_length, trailing_ref_pos = trailing_edge_poor_alignment_length_short_short(
            align, read_seq, ref_file, contiguous_match_count)
    return leading_length, leading_ref_pos, trailing_length, trailing_ref_pos


def get_sv_breakend_candidate_semi_aligned_short(
        read, align, ref_file, min_baseq=20, min_high_baseq_frac=0.8, contiguous_match_count=5):


    leading_edge_len = 0
    leading_edge_ref_pos = align["ref_start"]
    trailing_edge_len = 0
    trailing_edge_ref_pos = align["ref_end"]


    is_overlapping = is_overlapping_pair_short(read)

    read_seq = read.query_sequence
    if not read_seq:
        return leading_edge_len, leading_edge_ref_pos, trailing_edge_len, trailing_edge_ref_pos

    quals = read.query_qualities
    if not quals:
        quals = [30] * len(read_seq)

    read_size = len(read_seq)


    matched_align = matchify_edge_soft_clip_short(align)


    leading_temp, leading_ref_pos = 0, align["ref_start"]
    trailing_temp, trailing_ref_pos = 0, align["ref_end"]

    if align['softclip_left'] > 0:
        leading_temp, leading_ref_pos = leading_edge_poor_alignment_length_short_short(
            matched_align, read_seq, ref_file, contiguous_match_count)

    if align['softclip_right'] > 0:
        trailing_temp, trailing_ref_pos = trailing_edge_poor_alignment_length_short_short(
            matched_align, read_seq, ref_file, contiguous_match_count)


    if leading_temp + trailing_temp >= read_size:
        return leading_edge_len, leading_edge_ref_pos, trailing_edge_len, trailing_edge_ref_pos


    if leading_temp > 0:

        if (not is_overlapping or
                read.is_supplementary or
                not read.is_reverse):

            high_q_count = sum(1 for i in range(min(leading_temp, len(quals)))
                               if quals[i] >= min_baseq)
            if leading_temp > 0 and (high_q_count / leading_temp) >= min_high_baseq_frac:
                leading_edge_len = leading_temp
                leading_edge_ref_pos = leading_ref_pos

    # 检查尾端质量
    if trailing_temp > 0:
        if (not is_overlapping or
                read.is_supplementary or
                read.is_reverse):

            start_pos = max(0, len(quals) - trailing_temp)
            high_q_count = sum(1 for i in range(start_pos, len(quals))
                               if quals[i] >= min_baseq)
            if trailing_temp > 0 and (high_q_count / trailing_temp) >= min_high_baseq_frac:
                trailing_edge_len = trailing_temp
                trailing_edge_ref_pos = trailing_ref_pos
    return leading_edge_len, leading_edge_ref_pos, trailing_edge_len, trailing_edge_ref_pos


def get_sv_candidates_from_semi_aligned_short(opt, dopt, read, ref_chr_length, sample_name, align,
                                        ref_file, candidates):
    result = get_sv_breakend_candidate_semi_aligned_short(read, align, ref_file, opt["minBasecallQuality"],
                                                    opt["minHighBasecallQualityFraction"])

    if result is None:
        return  #
    leading_mismatch_len, leading_ref_pos, trailing_mismatch_len, trailing_ref_pos = result

    if (leading_mismatch_len + trailing_mismatch_len) >= len(read.query_sequence):
        return
    sv_evidence_source = SVET_SEMIALIGN
    is_complex = True
    before_breakend_offset = dopt.get("beforeBreakend", 100)
    after_breakend_offset = dopt.get("afterBreakend", 100)
    # 只输出有断点的一端
    sv_evidence_type = 'SEMI'
    if leading_mismatch_len >= opt["minSemiAlignedMismatchLen"]:
        pos = leading_ref_pos
        #print(f'leading_mismatch_len:{leading_mismatch_len}, leading_ref_pos:{leading_ref_pos}')
        candidates.append(get_split_sv_candidate_short(
            align, ref_chr_length, sample_name, pos, pos,
            before_breakend_offset, after_breakend_offset,
            sv_evidence_source, sv_evidence_type, is_complex))
    if trailing_mismatch_len >= opt["minSemiAlignedMismatchLen"]:
        pos = trailing_ref_pos
        candidates.append(get_split_sv_candidate_short(
            align, ref_chr_length, sample_name, pos, pos,
            before_breakend_offset, after_breakend_offset,
            sv_evidence_source, sv_evidence_type, is_complex))


def extract_semi_aligned_short(read, ref_file, ref_chr_length, sample_name, min_mismatch_len=30, bamfile=None):
    candidates = []
    if read.is_unmapped or not read.query_sequence:
        return candidates
    if read.has_tag('SA'):
        return candidates
    chr_name = read.reference_name
    ref_start = read.reference_start
    read_name = read.query_name
    cigarstring = read.cigarstring
    is_reverse = read.is_reverse
    primary_read = parse_cigar_string_short_save_string(chr_name, cigarstring, ref_start, read_name, read.mapping_quality,
                                                  is_reverse)
    dopt = {
        "beforeBreakend": 100,
        "afterBreakend": 100,
        "isSmallCandidates": True,
        "isTranscriptStrandKnown": False
    }
    is_read2 = read.is_paired and not read.is_read1
    frag_source = 'read2' if is_read2 else 'read1'
    opt = {
        "minSemiAlignedMismatchLen": min_mismatch_len,
        "useOverlapPairEvidence": True,
        "minBasecallQuality": 20,
        "minHighBasecallQualityFraction": 0.8
    }
    get_sv_candidates_from_semi_aligned_short(opt, dopt, read, ref_chr_length, sample_name, primary_read,
                                        ref_file, candidates)
    return candidates


def parse_cigar_string_short_save_string(chr_name, cigarstring, ref_start, read_name, mapping_quality, is_reverse=False):
    result = {
        'chr': chr_name,
        'softclip_left': 0,
        'softclip_right': 0,
        'match_length': 0,
        'insert_length': 0,
        'delete_length': 0,
        'mismatch_length': 0,
        'ref_start': ref_start,
        'ref_end': ref_start,
        'query_start': 0,
        'query_end': 0,
        'read_name': read_name,
        'is_reverse': is_reverse,
        'mapping_quality': mapping_quality,
        'cigarstring': None
    }
    cigar_ops = re.findall(r'(\d+)([MIDNSHP=X])', cigarstring)
    # For reverse aligned reads, CIGAR operations need to be reversed for parse, but NOT for path traversal!
    # if is_reverse:
    #     cigar_ops = cigar_ops[::-1]
    query_pos = 0
    ref_pos = ref_start
    for i, (length_str, op) in enumerate(cigar_ops):
        length = int(length_str)
        if op == 'M':
            result['match_length'] += length
            query_pos += length
            ref_pos += length
        elif op == 'I':
            result['insert_length'] += length
            query_pos += length
        elif op == 'D':
            result['delete_length'] += length
            ref_pos += length
        elif op == 'S':
            if i == 0:
                result['softclip_left'] = length
            elif i == len(cigar_ops) - 1:
                result['softclip_right'] = length
            query_pos += length
        elif op == '=':
            result['match_length'] += length
            query_pos += length
            ref_pos += length
        elif op == 'X':
            result['mismatch_length'] += length
            query_pos += length
            ref_pos += length
        elif op == 'N':
            # N is a skip operation, it does not affect query or reference positions
            ref_pos += length
        elif op == 'H':
            pass  # Hard clipping does not affect query or reference positions
    result['ref_end'] = ref_pos
    result['query_start'] = result['softclip_left']
    result['query_end'] = query_pos - result['softclip_right']
    result['cigarstring'] = cigarstring_to_path_short(cigarstring)
    return result


def save_contig_to_fastq(contig_seq, output_fastq, contig_id="assembled_contig"):

    with open(output_fastq, "a") as fq:
        fq.write(f"@{contig_id}\n")
        fq.write(f"{contig_seq}\n")
        fq.write("+\n")
        fq.write("I" * len(contig_seq) + "\n")  # 假设质量值为 'I'


def has_assemble_in_region(candidates, count, save_path, sv_type):

    if not candidates:
        print(f"No reads for assembly in {sv_type}")
        return False, "No reads for assembly"

    try:

        min_k, max_k = 35, 75
        contig_num = 10
        assembled_read_lis = mantaAss(candidates, min_k, max_k, contig_num)


        contig_count = 0
        for i, contig_seq in enumerate(assembled_read_lis[1]):
            if not contig_seq:
                continue
            contig_name = f'contig_{count}_{sv_type}_{i}'

            if isinstance(contig_seq, bytes):
                contig_seq_str = contig_seq.decode()
            else:
                contig_seq_str = str(contig_seq)

            save_contig_to_fastq(contig_seq_str, save_path, contig_id=contig_name)
            contig_count += 1

        print(f"Assembled {contig_count} contigs for {sv_type}")
        return True, f"Assembled {contig_count} contigs"

    except Exception as e:
        print(f"Assembly failed for {sv_type}: {e}")
        return False, f"Assembly failed: {e}"

def collect_evidence_reads_small_short(bam_normal, chrom_start, start, end, min_mapq, ref_chr_length,
                                       sample_name, chrom_list, ref_chr_id_to_name, ref_path,
                                       min_indel=40, max_isize=1000000, breakpoint_margin=500):

    import pysam
    candidates = []
    start = start - breakpoint_margin
    end = end + breakpoint_margin
    bam_normal_mate = pysam.AlignmentFile(bam_normal.filename, "rb")

    try:
        for read in bam_normal.fetch(chrom_start, start, end):
            if is_read_filtered_short(read, min_mapq):
                continue

            if read.reference_end is None:
                continue

            if not (read.reference_start < end and read.reference_end > start):
                continue

            one = False

            # Split reads
            split_candidates = extract_split_reads_short(read, ref_chr_length, sample_name, chrom_list, min_mapq)

            # Large indels
            indel_candidates = get_sv_candidates_from_read_indels_short(read, ref_chr_length, sample_name, start, end, min_size=min_indel)

            # Discordant pairs
            discordant_candidates = extract_discordant_pairs_short(
                read, ref_chr_length, sample_name, ref_chr_id_to_name,
                min_isize=500, max_isize=max_isize
            )

            # Semi-aligned
            semi_candidates = []
            if read.reference_end and read.reference_start < read.reference_end:
                semi_candidates = extract_semi_aligned_short(read, ref_path, ref_chr_length, sample_name, min_mismatch_len=20)

            if split_candidates and not one:
                candidates.append(read.query_sequence.encode())
                one = True
            elif discordant_candidates and not one:
                candidates.append(read.query_sequence.encode())
                try:
                    if read.is_paired and not read.mate_is_unmapped:
                        mate_chrom = bam_normal_mate.get_reference_name(read.next_reference_id)
                        mate_pos = read.next_reference_start
                        for mate_read in bam_normal_mate.fetch(mate_chrom, max(0, mate_pos - 100), mate_pos + 100):
                            if mate_read.query_name == read.query_name:
                                candidates.append(mate_read.query_sequence.encode())
                                break
                except Exception as e:
                    print(f"Failed to find mate: {e}")
                one = True
            elif indel_candidates and not one:
                candidates.append(read.query_sequence.encode())
                one = True
            elif semi_candidates and not one:
                candidates.append(read.query_sequence.encode())
                one = True

    finally:
        bam_normal_mate.close()

    return candidates


def collect_evidence_reads_large_short(bam_normal, chrom, start, end, min_mapq, ref_chr_length,
                                       sample_name, chrom_list, ref_chr_id_to_name, ref_path,
                                       min_isize_large, max_isize_large, breakpoint_margin=500):

    import pysam
    candidates_start = []
    candidates_end = []

    bam_normal_mate = pysam.AlignmentFile(bam_normal.filename, "rb")
    try:
        for read in bam_normal.fetch(chrom, max(0, start - breakpoint_margin), start + breakpoint_margin):
            if is_read_filtered_short(read, min_mapq):
                continue

            one = False

            split_candidates = extract_split_reads_short(read, ref_chr_length, sample_name, chrom_list, min_mapq)

            if read.cigartuples and any(op == 4 and length >= 15 for op, length in read.cigartuples):
                if not one:
                    candidates_start.append(read.query_sequence.encode())
                    one = True

            discordant_candidates = extract_discordant_pairs_short(
                read, ref_chr_length, sample_name, ref_chr_id_to_name,
                min_isize=min_isize_large, max_isize=max_isize_large
            )

            if split_candidates and not one:
                candidates_start.append(read.query_sequence.encode())
                one = True
            elif discordant_candidates and not one:
                candidates_start.append(read.query_sequence.encode())

                try:
                    if read.is_paired and not read.mate_is_unmapped:
                        mate_chrom = bam_normal_mate.get_reference_name(read.next_reference_id)
                        mate_pos = read.next_reference_start
                        for mate_read in bam_normal_mate.fetch(mate_chrom, max(0, mate_pos - 100), mate_pos + 100):
                            if mate_read.query_name == read.query_name:
                                candidates_start.append(mate_read.query_sequence.encode())
                                break
                except Exception as e:
                    print(f"Failed to find mate: {e}")
                one = True

        for read in bam_normal.fetch(chrom, max(0, end - breakpoint_margin), end + breakpoint_margin):
            if is_read_filtered_short(read, min_mapq):
                continue

            one = False
            split_candidates = extract_split_reads_short(read, ref_chr_length, sample_name, chrom_list, min_mapq)

            if read.cigartuples and any(op == 4 and length >= 15 for op, length in read.cigartuples):
                if not one:
                    candidates_end.append(read.query_sequence.encode())
                    one = True

            discordant_candidates = extract_discordant_pairs_short(
                read, ref_chr_length, sample_name, ref_chr_id_to_name,
                min_isize=min_isize_large, max_isize=max_isize_large
            )

            if split_candidates and not one:
                candidates_end.append(read.query_sequence.encode())
                one = True
            elif discordant_candidates and not one:
                candidates_end.append(read.query_sequence.encode())

                try:
                    if read.is_paired and not read.mate_is_unmapped:
                        mate_chrom = bam_normal_mate.get_reference_name(read.next_reference_id)
                        mate_pos = read.next_reference_start
                        for mate_read in bam_normal_mate.fetch(mate_chrom, max(0, mate_pos - 100), mate_pos + 100):
                            if mate_read.query_name == read.query_name:
                                candidates_end.append(mate_read.query_sequence.encode())
                                break
                except Exception as e:
                    print(f"Failed to find mate: {e}")
                one = True

    finally:
        bam_normal_mate.close()

    return candidates_start, candidates_end


def collect_evidence_reads_BND_short(bam_normal, chrom_start, chrom_end, start, end, min_mapq,
                                     ref_chr_length, sample_name, chrom_list, ref_chr_id_to_name,
                                     ref_path, min_isize_large, max_isize_large, breakpoint_margin=500):

    import pysam
    candidates_start = []
    candidates_end = []
    bam_normal_mate = pysam.AlignmentFile(bam_normal.filename, "rb")

    try:
        for read in bam_normal.fetch(chrom_start, max(0, start - breakpoint_margin), start + breakpoint_margin):
            if is_read_filtered_short(read, min_mapq):
                continue

            one = False
            split_candidates = extract_split_reads_short(read, ref_chr_length, sample_name, chrom_list, min_mapq)

            if read.cigartuples and any(op == 4 and length >= 15 for op, length in read.cigartuples):
                if not one:
                    candidates_start.append(read.query_sequence.encode())
                    one = True

            if (read.is_paired and not read.mate_is_unmapped and
                read.reference_id != read.next_reference_id):
                if not one:
                    candidates_start.append(read.query_sequence.encode())
                    try:
                        mate_chrom = bam_normal_mate.get_reference_name(read.next_reference_id)
                        mate_pos = read.next_reference_start
                        for mate_read in bam_normal_mate.fetch(mate_chrom, max(0, mate_pos - 100), mate_pos + 100):
                            if mate_read.query_name == read.query_name:
                                candidates_start.append(mate_read.query_sequence.encode())
                                break
                    except Exception as e:
                        print(f"Failed to find mate: {e}")
                    one = True

            if split_candidates and not one:
                candidates_start.append(read.query_sequence.encode())
                one = True

        if chrom_start != chrom_end:
            for read in bam_normal.fetch(chrom_end, max(0, end - breakpoint_margin), end + breakpoint_margin):
                if is_read_filtered_short(read, min_mapq):
                    continue

                one = False
                split_candidates = extract_split_reads_short(read, ref_chr_length, sample_name, chrom_list, min_mapq)

                if read.cigartuples and any(op == 4 and length >= 15 for op, length in read.cigartuples):
                    if not one:
                        candidates_end.append(read.query_sequence.encode())
                        one = True

                if (read.is_paired and not read.mate_is_unmapped and
                    read.reference_id != read.next_reference_id):
                    if not one:
                        candidates_end.append(read.query_sequence.encode())
                        try:
                            mate_chrom = bam_normal_mate.get_reference_name(read.next_reference_id)
                            mate_pos = read.next_reference_start
                            for mate_read in bam_normal_mate.fetch(mate_chrom, max(0, mate_pos - 100), mate_pos + 100):
                                if mate_read.query_name == read.query_name:
                                    candidates_end.append(mate_read.query_sequence.encode())
                                    break
                        except Exception as e:
                            print(f"Failed to find mate: {e}")
                        one = True

                if split_candidates and not one:
                    candidates_end.append(read.query_sequence.encode())
                    one = True

    finally:
        bam_normal_mate.close()

    return candidates_start, candidates_end
