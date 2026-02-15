import re
import numpy as np
from numba import jit


def parse_cigar_string(chr_name, cigarstring, ref_start, read_name, mapping_quality, is_reverse=False):
    """
    Parse CIGAR string, calculate reference genome alignment start, end position and query sequence position.
    """
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

    cigar_ops = re.findall(r'(\d+)([MIDNSHP=X])', cigarstring)

    if is_reverse:
        cigar_ops = cigar_ops[::-1]

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
            else:
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

    result['ref_end'] = ref_pos
    result['query_start'] = result['softclip_left']
    result['query_end'] = query_pos - result['softclip_right']

    return result


@jit(nopython=True)
def merge_intervals_numba(intervals, max_gap):
    if len(intervals) == 0:
        return np.empty((0, 3), dtype=np.int64)

    sorted_indices = np.argsort(intervals[:, 0])
    sorted_intervals = intervals[sorted_indices]

    merged = []
    current_start = sorted_intervals[0, 0]
    current_end = sorted_intervals[0, 1]
    current_length = sorted_intervals[0, 2]

    for i in range(1, len(sorted_intervals)):
        start = sorted_intervals[i, 0]
        end = sorted_intervals[i, 1]
        length = sorted_intervals[i, 2]

        if start - current_end <= max_gap:
            current_end = max(current_end, end)
            current_length += length
        else:
            merged.append([current_start, current_end, current_length])
            current_start, current_end, current_length = start, end, length

    merged.append([current_start, current_end, current_length])
    return np.array(merged)


@jit(nopython=True)
def process_cigar_numba(cigar_ops, ref_start, min_size, query_sequence_length):
    del_intervals = []
    ins_intervals = []

    reference_pos = ref_start
    query_pos = 0

    for i in range(len(cigar_ops)):
        operation = cigar_ops[i, 0]
        length = cigar_ops[i, 1]

        if operation in (0, 7, 8):
            reference_pos += length
            query_pos += length
        elif operation == 2:
            if length >= min_size:
                del_intervals.append([reference_pos, reference_pos + length, length])
            reference_pos += length
        elif operation == 1:
            if length >= min_size:
                ins_intervals.append([reference_pos, reference_pos, length, query_pos, query_pos + length])
            query_pos += length

    return del_intervals, ins_intervals


def process_cigar_operations_long(read, start, end, platform='HIFI', min_size=30, max_gap=30):
    """
    Process the CIGAR string from a read, identifying large insertions and deletions.
    Returns two lists: deletions and insertions, where insertions include the actual sequence.
    """
    del_intervals = []
    ins_intervals = []

    cigar_array = np.array(read.cigartuples, dtype=np.int32)
    query_len = (
        len(read.query_sequence)
        if read.query_sequence
        else (read.query_length if read.query_length else 0)
    )
    del_numba, ins_numba = process_cigar_numba(cigar_array, read.reference_start, min_size, query_len)

    for del_interval in del_numba:
        ref_start, ref_end, length = del_interval
        if start <= ref_end and end >= ref_start:
            del_intervals.append([ref_start, ref_end, length])

    for ins_interval in ins_numba:
        ref_pos, _, length, query_start, query_end = ins_interval
        if start <= ref_pos and end >= ref_pos:
            insert_seq = (read.query_sequence[query_start:query_end]
                          if read.query_sequence else None)
            ins_intervals.append([ref_pos, ref_pos, length, insert_seq])
    if platform == 'ONT':
        del_intervals = merge_intervals(del_intervals, max_gap=max_gap)
        ins_intervals = merge_intervals(ins_intervals, max_gap=max_gap, is_insertion=True)
    else:
        del_intervals = merge_intervals(del_intervals, max_gap=max_gap)
        ins_intervals = merge_intervals(ins_intervals, max_gap=max_gap, is_insertion=True)

    return del_intervals, ins_intervals


def merge_intervals(intervals, max_gap=0, is_insertion=False):
    if not intervals:
        return []

    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    merged = [sorted_intervals[0]]

    for current in sorted_intervals[1:]:
        last = merged[-1]

        if current[0] <= last[1] + max_gap:
            merged[-1] = [
                last[0],
                max(last[1], current[1]),
                last[2] + current[2]
            ]
            if is_insertion and len(current) > 3:
                if len(merged[-1]) <= 3:
                    merged[-1].append(current[3])
                else:
                    merged[-1][3] += current[3]
        else:
            merged.append(current)

    return merged


def process_cigar_operations_refine(read, start, end, min_size=30, max_gap=150):
    """
    Process the CIGAR string from a read, identifying large insertions and deletions.
    Returns two lists: deletions and insertions, where insertions include the actual sequence.
    """
    del_intervals = []
    ins_intervals = []
    reference_pos = read.reference_start
    query_pos = 0

    for operation, length in read.cigartuples:
        if operation in [0, 7, 8]:
            reference_pos += length
            query_pos += length
        elif operation == 2:
            if length >= min_size and start <= reference_pos + length and end >= reference_pos:
                del_intervals.append([reference_pos, reference_pos + length, length])
            reference_pos += length
        elif operation == 1:
            if length >= min_size and start <= reference_pos and end >= reference_pos:
                if read.query_sequence is None:
                    insert_seq = None
                else:
                    insert_seq = read.query_sequence[query_pos:query_pos + length]
                ins_intervals.append([reference_pos, reference_pos, length, insert_seq])
            query_pos += length

    del_intervals = merge_intervals(del_intervals, max_gap=max_gap)
    ins_intervals = merge_intervals(ins_intervals, max_gap=max_gap, is_insertion=True)

    return del_intervals, ins_intervals


def process_cigar_numba1(cigar_ops, ref_start, min_size, query_sequence_length):
    del_intervals = []
    ins_intervals = []

    reference_pos = ref_start
    query_pos = 0

    for i in range(len(cigar_ops)):
        operation = cigar_ops[i, 1]
        length = cigar_ops[i, 0]

        if operation in (0, 7, 8):
            reference_pos += length
            query_pos += length
        elif operation == 2:
            if length >= min_size:
                del_intervals.append([reference_pos, reference_pos + length, length])
            reference_pos += length
        elif operation == 1:
            if length >= min_size:
                ins_intervals.append([reference_pos, reference_pos, length, query_pos, query_pos + length])
            query_pos += length

    return del_intervals, ins_intervals