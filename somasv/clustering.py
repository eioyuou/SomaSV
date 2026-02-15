import numpy as np
import pysam
import multiprocessing
from numba import jit, njit
from collections import defaultdict

from somasv.utils import calculate_statistics_numba


# ==================== Cluster Initialization and Management ====================

@njit
def check_overlap_numba(cluster_start, cluster_end, breakpoint_pos, extension):
    distance = breakpoint_pos - cluster_end
    return cluster_start <= breakpoint_pos and distance <= extension


def is_overlap_with_cluster(cluster, breakpoint, extension):
    cluster_start = cluster['start']
    cluster_end = cluster['end']
    breakpoint_pos = breakpoint['start_loc']
    return check_overlap_numba(cluster_start, cluster_end, breakpoint_pos, extension)


@njit
def update_cluster_positions_numba(cluster_start, cluster_end, new_start, new_end, is_nearby_event):
    if cluster_start <= cluster_end:
        updated_start = min(cluster_start, new_start)
        if is_nearby_event:
            updated_end = max(cluster_end, new_end)
        else:
            updated_end = max(cluster_end, new_start)
    else:
        updated_start = max(cluster_start, new_start)
        if is_nearby_event:
            updated_end = min(cluster_end, new_end)
        else:
            updated_end = min(cluster_end, new_start)

    return updated_start, updated_end


def initialize_cluster(initial_breakpoint):
    """ Initialize a new breakpoint cluster dictionary """
    cluster = {
        "chr": initial_breakpoint['start_chr'],
        "start": initial_breakpoint['start_loc'],
        "end": initial_breakpoint['end_loc'] if initial_breakpoint['is_nearby_event'] else initial_breakpoint[
            'start_loc'],
        "source": initial_breakpoint['source'],
        "breakpoints": [initial_breakpoint],
        "supporting_reads": {initial_breakpoint['read_name']},
        "stats": None
    }
    return cluster


def merge_breakpoint_into_cluster(cluster, new_breakpoint):
    new_start, new_end = update_cluster_positions_numba(
        cluster['start'],
        cluster['end'],
        new_breakpoint['start_loc'],
        new_breakpoint['end_loc'] if new_breakpoint['is_nearby_event'] else new_breakpoint['start_loc'],
        new_breakpoint['is_nearby_event']
    )

    cluster['start'] = new_start
    cluster['end'] = new_end
    cluster['breakpoints'].append(new_breakpoint)

    if new_breakpoint['read_name'] not in cluster['supporting_reads']:
        cluster['supporting_reads'].add(new_breakpoint['read_name'])


def calculate_cluster_statistics(cluster):
    if not cluster['stats']:
        starts = []
        mapqs = []
        event_sizes = []

        for bp in cluster['breakpoints']:
            starts.append(bp['start_loc'])
            mapqs.append(bp['mapq'])
            if bp['breakpoint_notation'] in ["<INS>", "+", "-"]:
                event_sizes.append(bp['insert_length'])
            elif bp['start_chr'] == bp['end_chr']:
                event_sizes.append(abs(bp['start_loc'] - bp['end_loc']))
            else:
                event_sizes.append(0)

        starts_std, mapq_mean, event_size_std, event_size_median, event_size_mean = \
            calculate_statistics_numba(starts, mapqs, event_sizes)

        cluster['stats'] = {
            'starts_std_dev': round(starts_std, 3),
            'mapq_mean': round(mapq_mean, 3),
            'event_size_std_dev': round(event_size_std, 3),
            'event_size_median': round(event_size_median, 3),
            'event_size_mean': round(event_size_mean, 3)
        }

    return cluster['stats']


def recalculate_supporting_reads(cluster):
    cluster['supporting_reads'] = {bp['read_name'] for bp in cluster['breakpoints']}
    return cluster


def cluster_by_insert_length(breakpoints, fraction):
    """
    Cluster breakpoints by insertion length.
    """
    insertion_breakpoints, single_breakpoints = [], []

    for bp in breakpoints:
        if bp['breakpoint_notation'] == "<INS>":
            insertion_breakpoints.append(bp)
        else:
            single_breakpoints.append(bp)

    insertion_breakpoints.sort(key=lambda bp: bp['insert_length'])

    from math import ceil
    clusters = []

    if insertion_breakpoints:
        for bp in insertion_breakpoints:
            if not clusters:
                clusters.append([bp])
            else:
                last_insert_length = clusters[-1][-1]['insert_length']
                threshold = ceil(last_insert_length * (2 - fraction))

                if bp['insert_length'] <= threshold:
                    clusters[-1].append(bp)
                else:
                    clusters.append([bp])

        clusters[-1].extend(single_breakpoints)
    else:
        clusters.append(single_breakpoints)

    return clusters


# ==================== End Cluster Functions ====================

def merge_end_breakpoint_into_cluster(cluster, new_breakpoint):
    """
    Add new breakpoint to cluster based on end position.
    """
    if cluster['start'] < cluster['end']:
        cluster['end'] = max(cluster['end'], new_breakpoint['end_loc'])
    else:
        cluster['end'] = min(cluster['end'], new_breakpoint['end_loc'])

    breakpoint_start = new_breakpoint['start_loc'] if new_breakpoint['is_nearby_event'] else new_breakpoint['end_loc']

    if cluster['start'] < cluster['end']:
        cluster['start'] = min(cluster['start'], breakpoint_start)
    else:
        cluster['start'] = max(cluster['start'], breakpoint_start)

    cluster['breakpoints'].append(new_breakpoint)
    cluster['supporting_reads'].add(new_breakpoint['read_name'])

    cluster['stats'] = None


def is_breakpoint_in_end_cluster(cluster, breakpoint, buffer):
    """ Check if breakpoint can join existing cluster (based on end position). """
    cluster_extended_start = cluster['start'] - buffer
    return cluster_extended_start <= breakpoint['end_loc']


def initialize_cluster_end(initial_breakpoint):
    """ Initialize a new breakpoint cluster dictionary """
    cluster = {
        "chr": initial_breakpoint['end_chr'],
        "start": initial_breakpoint['end_loc'],
        "end": initial_breakpoint['start_loc'] if initial_breakpoint['is_nearby_event'] else initial_breakpoint[
            'end_loc'],
        "source": initial_breakpoint['source'],
        "breakpoints": [initial_breakpoint],
        "supporting_reads": {initial_breakpoint['read_name']},
        "stats": None
    }
    return cluster


# ==================== CIGAR Refinement (integrated) ====================

def has_significant_variations(cigar_info):
    if not cigar_info:
        return False

    significant_vars = 0
    for del_interval in cigar_info['deletions']:
        if del_interval[2] >= 30:
            significant_vars += 1

    for ins_interval in cigar_info['insertions']:
        if ins_interval[2] >= 30:
            significant_vars += 1

    return significant_vars > 0


@njit
def calculate_interval_overlap_optimized_numba(intervals1, intervals2, tolerance):
    overlap_count = 0
    for i in range(intervals1.shape[0]):
        start1, end1, len1 = intervals1[i, 0], intervals1[i, 1], intervals1[i, 2]
        for j in range(intervals2.shape[0]):
            start2, end2, len2 = intervals2[j, 0], intervals2[j, 1], intervals2[j, 2]
            start_diff = abs(start1 - start2)
            end_diff = abs(end1 - end2)
            if start_diff <= tolerance and end_diff <= tolerance:
                min_len = min(len1, len2)
                max_len = max(len1, len2)
                if max_len > 0:
                    size_ratio = min_len / max_len
                    if size_ratio >= 0.8:
                        overlap_count += 1
                        break
    return overlap_count


def calculate_interval_overlap_simple(intervals1, intervals2, tolerance):
    overlap_count = 0
    for iv1 in intervals1:
        if len(iv1) < 3:
            continue
        start1, end1, len1 = iv1[0], iv1[1], iv1[2]
        for iv2 in intervals2:
            if len(iv2) < 3:
                continue
            start2, end2, len2 = iv2[0], iv2[1], iv2[2]
            if (abs(start1 - start2) <= tolerance and
                    abs(end1 - end2) <= tolerance and
                    min(len1, len2) / max(len1, len2) >= 0.8):
                overlap_count += 1
                break
    return overlap_count


def calculate_interval_overlap_optimized(intervals1, intervals2, tolerance):
    if not intervals1 or not intervals2:
        return 0
    try:
        arr1 = np.array([[iv[0], iv[1], iv[2]] for iv in intervals1], dtype=np.float64)
        arr2 = np.array([[iv[0], iv[1], iv[2]] for iv in intervals2], dtype=np.float64)
        return calculate_interval_overlap_optimized_numba(arr1, arr2, tolerance)
    except (IndexError, ValueError):
        return calculate_interval_overlap_simple(intervals1, intervals2, tolerance)


def calculate_cigar_similarity_optimized(cigar1, cigar2, position_tolerance=50):
    if not cigar1 or not cigar2:
        return 0.0

    del1, ins1 = cigar1.get('deletions', []), cigar1.get('insertions', [])
    del2, ins2 = cigar2.get('deletions', []), cigar2.get('insertions', [])

    if not (del1 or ins1) and not (del2 or ins2):
        return 1.0
    if not (del1 or ins1) or not (del2 or ins2):
        return 0.0

    del_overlap = calculate_interval_overlap_optimized(del1, del2, position_tolerance)
    ins_overlap = calculate_interval_overlap_optimized(ins1, ins2, position_tolerance)

    total_vars1 = len(del1) + len(ins1)
    total_vars2 = len(del2) + len(ins2)

    if total_vars1 == 0 or total_vars2 == 0:
        return 0.0

    total_overlap = del_overlap + ins_overlap
    avg_total_vars = (total_vars1 + total_vars2) / 2.0

    similarity = total_overlap / avg_total_vars
    return min(similarity, 1.0)


@njit
def calculate_overlap_score_numba(intervals1, intervals2, tolerance):
    overlap_count = 0
    for i in range(intervals1.shape[0]):
        for j in range(intervals2.shape[0]):
            pos_diff1 = abs(intervals1[i, 0] - intervals2[j, 0])
            pos_diff2 = abs(intervals1[i, 1] - intervals2[j, 1])
            if pos_diff1 <= tolerance and pos_diff2 <= tolerance:
                size1, size2 = intervals1[i, 2], intervals2[j, 2]
                size_ratio = min(size1, size2) / max(size1, size2) if max(size1, size2) > 0 else 0
                if size_ratio >= 0.8:
                    overlap_count += 1
                    break
    return overlap_count


def calculate_cigar_similarity(cigar1, cigar2, position_tolerance=50):
    if not cigar1 or not cigar2:
        return 0.0

    def to_numpy_intervals(intervals):
        if not intervals:
            return np.empty((0, 3), dtype=np.float64)
        return np.array([[iv[0], iv[1], iv[2]] for iv in intervals], dtype=np.float64)

    del1_array = to_numpy_intervals(cigar1['deletions'])
    del2_array = to_numpy_intervals(cigar2['deletions'])
    ins1_array = to_numpy_intervals(cigar1['insertions'])
    ins2_array = to_numpy_intervals(cigar2['insertions'])

    del_overlap = calculate_overlap_score_numba(del1_array, del2_array, position_tolerance)
    ins_overlap = calculate_overlap_score_numba(ins1_array, ins2_array, position_tolerance)

    total_vars1 = len(cigar1['deletions']) + len(cigar1['insertions'])
    total_vars2 = len(cigar2['deletions']) + len(cigar2['insertions'])

    if total_vars1 == 0 and total_vars2 == 0:
        return 1.0
    if total_vars1 == 0 or total_vars2 == 0:
        return 0.0

    overlap_count = del_overlap + ins_overlap
    avg_total_vars = (total_vars1 + total_vars2) / 2

    return overlap_count / avg_total_vars


@njit
def calculate_interval_overlap_numba(intervals1_array, intervals2_array, tolerance):
    overlap_count = 0
    for i in range(intervals1_array.shape[0]):
        for j in range(intervals2_array.shape[0]):
            if (abs(intervals1_array[i, 0] - intervals2_array[j, 0]) <= tolerance and
                    abs(intervals1_array[i, 1] - intervals2_array[j, 1]) <= tolerance):
                size_ratio = min(intervals1_array[i, 2], intervals2_array[j, 2]) / max(intervals1_array[i, 2], intervals2_array[j, 2])
                if size_ratio >= 0.8:
                    overlap_count += 1
                    break
    return overlap_count


def calculate_interval_overlap(intervals1, intervals2, tolerance):
    if not intervals1 or not intervals2:
        return 0
    intervals1_array = np.array([[int1[0], int1[1], int1[2]] for int1 in intervals1], dtype=np.float64)
    intervals2_array = np.array([[int2[0], int2[1], int2[2]] for int2 in intervals2], dtype=np.float64)
    return calculate_interval_overlap_numba(intervals1_array, intervals2_array, tolerance)


def extract_cigar_info(breakpoints, bamfile_path, cluster, min_size=2):
    """Batch-optimized extraction of CIGAR information"""
    cigar_info = {}
    read_name_to_bp = {bp['read_name']: bp for bp in breakpoints}
    target_reads = set(read_name_to_bp.keys())

    chr_name = cluster['chr']
    start_pos = max(0, cluster['start'] - 500)
    end_pos = cluster['end'] + 500

    with pysam.AlignmentFile(bamfile_path, "rb") as bamfile:
        reads_batch = []
        for read in bamfile.fetch(chr_name, start_pos, end_pos):
            if read.query_name in target_reads:
                reads_batch.append(read)

        cigar_info = process_cigars_batch(reads_batch, read_name_to_bp, min_size)

    return cigar_info


def process_cigars_batch(reads_batch, read_name_to_bp, min_size):
    """Batch processing of CIGAR strings"""
    cigar_info = {}

    for read in reads_batch:
        bp = read_name_to_bp[read.query_name]

        del_intervals, ins_intervals = process_cigar_operations_optimized(
            read, min_size
        )

        cigar_info[read.query_name] = {
            'breakpoint': bp,
            'deletions': del_intervals,
            'insertions': ins_intervals,
            'cigar_string': read.cigarstring,
            'reference_start': read.reference_start,
            'reference_end': read.reference_end
        }

    return cigar_info


def union_find_clustering(similarity_matrix, threshold, n_reads):
    parent = list(range(n_reads))

    def find(x):
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    for (i, j), similarity in similarity_matrix.items():
        if similarity >= threshold:
            union(i, j)

    clusters = defaultdict(list)
    for i in range(n_reads):
        root = find(i)
        clusters[root].append(i)

    return list(clusters.values())


def convert_to_bp_clusters(index_clusters, read_names, breakpoints):
    read_name_to_bp = {bp['read_name']: bp for bp in breakpoints}

    bp_clusters = []
    for cluster_indices in index_clusters:
        cluster_bps = []
        for idx in cluster_indices:
            read_name = read_names[idx]
            if read_name in read_name_to_bp:
                cluster_bps.append(read_name_to_bp[read_name])
        if cluster_bps:
            bp_clusters.append(cluster_bps)

    return bp_clusters


def cluster_by_cigar_similarity(cigar_info, breakpoints, threshold):
    if not cigar_info:
        return [breakpoints]

    read_names = list(cigar_info.keys())
    n_reads = len(read_names)

    if n_reads < 2:
        return [breakpoints]

    similarity_matrix = {}
    for i in range(n_reads):
        for j in range(i + 1, n_reads):
            read1, read2 = read_names[i], read_names[j]
            similarity = calculate_cigar_similarity_optimized(
                cigar_info[read1], cigar_info[read2]
            )
            similarity_matrix[(i, j)] = similarity

    clusters = union_find_clustering(similarity_matrix, threshold, n_reads)

    return convert_to_bp_clusters(clusters, read_names, breakpoints)


def create_refined_cluster_simple(original_cluster, all_breakpoints):
    """ Create a simplified refined cluster """
    if len(all_breakpoints) >= 1:
        refined_cluster = {
            'chr': original_cluster['chr'],
            'start': min(bp['start_loc'] for bp in all_breakpoints),
            'end': max(bp['end_loc'] if bp.get('end_loc') else bp['start_loc'] for bp in all_breakpoints),
            'source': original_cluster['source'],
            'breakpoints': all_breakpoints,
            'supporting_reads': {bp['read_name'] for bp in all_breakpoints},
            'stats': None
        }
        return refined_cluster
    return None


def refine_cluster_by_cigar(cluster, tumor_bamfile, normal_bamfile, similarity_threshold=0.6):
    """
    Refine a cluster based on CIGAR strings.
    """
    normal_bp = [bp for bp in cluster['breakpoints'] if bp['sample'] == 'normal']
    tumor_bp = [bp for bp in cluster['breakpoints'] if bp['sample'] == 'tumor']

    if not tumor_bp:
        return [cluster]

    cluster_span = cluster['end'] - cluster['start']
    if cluster_span < 100:
        return [cluster]

    tumor_cigar_info = extract_cigar_info(tumor_bp, tumor_bamfile, cluster)

    valid_tumor_cigar = {k: v for k, v in tumor_cigar_info.items()
                         if has_significant_variations(v)}

    if len(valid_tumor_cigar) < 2:
        return [cluster]

    tumor_refined_clusters = cluster_by_cigar_similarity(tumor_cigar_info, tumor_bp, similarity_threshold)

    refined_clusters = []
    for tumor_cluster_bps in tumor_refined_clusters:
        all_breakpoints = tumor_cluster_bps + normal_bp
        refined_cluster = create_refined_cluster_simple(original_cluster=cluster,
                                                        all_breakpoints=all_breakpoints)
        if refined_cluster:
            refined_clusters.append(refined_cluster)

    return refined_clusters if refined_clusters else [cluster]


# ==================== Group Related Breakpoints ====================

def group_related_breakpoints(chromosome, breakpoints, tumor_bamfile, normal_bamfile, extension,
                              insertion_additional=None):
    """ Given a list of breakpoints, cluster them based on location and type """
    stack = []
    breakpoints.sort(key=lambda x: x['start_loc'])
    for bp in breakpoints:
        bp_notation_type = str(bp['breakpoint_notation'])
        bp_extension = insertion_additional if (insertion_additional and bp_notation_type == "<INS>") else extension

        if len(stack) == 0:
            new_cluster = initialize_cluster(bp)
            stack.append(new_cluster)
        elif not is_overlap_with_cluster(stack[-1], bp, bp_extension):
            new_cluster = initialize_cluster(bp)
            stack.append(new_cluster)
        else:
            merge_breakpoint_into_cluster(stack[-1], bp)

    filtered_stack = [cluster for cluster in stack if len(cluster['supporting_reads']) >= 2]

    for cluster in filtered_stack:
        calculate_cluster_statistics(cluster)
    print(f"Chromosome {chromosome} has {len(filtered_stack)} clusters after refinement.")
    return chromosome, filtered_stack


def group_related_breakpoints_task(args):
    """Multiprocess task function for clustering breakpoints"""
    print(f'Processing {args[0]}')
    chrom, breakpoints, tumor_bamfile, normal_bamfile, extension, insertion_additional = args
    return group_related_breakpoints(chrom, breakpoints, tumor_bamfile, normal_bamfile, extension,
                                     insertion_additional)


def parallel_group_related_breakpoints(breakpoints_by_chrom, tumor_bamfile, normal_bamfile, extension,
                                       insertion_additional=None, num_processes=None):
    """Use multiprocessing to cluster breakpoints for multiple chromosomes"""
    tasks = [(chrom, breakpoints, tumor_bamfile, normal_bamfile, extension, insertion_additional) for
             chrom, breakpoints in
             breakpoints_by_chrom.items()]

    with multiprocessing.Pool(processes=num_processes) as pool:
        results = pool.map(group_related_breakpoints_task, tasks)

    clustered_breakpoints = {chrom: clusters for chrom, clusters in results}

    return clustered_breakpoints


def group_related_end_breakpoints(chromosome, breakpoints, tumor_bamfile, normal_bamfile, extension,
                                  insertion_additional=None):
    """
    Cluster breakpoints based on end position (end_loc)
    """
    stack = []
    breakpoints.sort(key=lambda x: x['end_loc'])
    for bp in breakpoints:
        bp_notation_type = str(bp['breakpoint_notation'])
        bp_extension = insertion_additional if (insertion_additional and bp_notation_type == "<INS>") else extension

        if len(stack) == 0:
            new_cluster = initialize_cluster_end(bp)
            stack.append(new_cluster)
        elif not is_breakpoint_in_end_cluster(stack[-1], bp, bp_extension):
            new_cluster = initialize_cluster_end(bp)
            stack.append(new_cluster)
        else:
            merge_end_breakpoint_into_cluster(stack[-1], bp)

    filtered_stack = [cluster for cluster in stack if len(cluster['supporting_reads']) >= 2]

    refined_stack = []
    for cluster in filtered_stack:
        refined_clusters = refine_cluster_by_cigar(cluster, tumor_bamfile, normal_bamfile, similarity_threshold=0.7)
        refined_stack.extend(refined_clusters)

    for cluster in refined_stack:
        calculate_cluster_statistics(cluster)

    return chromosome, refined_stack