import numpy as np
from numba import jit
from math import ceil


@jit(nopython=True)
def add_coverage_numba(coverage_array, starts, ends, coverage_binsize):
    """numba加速的覆盖度累加"""
    for i in range(len(starts)):
        start_bin = starts[i] // coverage_binsize
        end_bin = ends[i] // coverage_binsize

        if start_bin > end_bin:
            start_bin, end_bin = end_bin, start_bin

        for bin_idx in range(start_bin, end_bin + 1):
            if 0 <= bin_idx < len(coverage_array):
                coverage_array[bin_idx] += 1


def add_coverage_to_arrays(coverage_arrays, start_positions, end_positions, haplotypes, coverage_binsize):
    if not isinstance(start_positions, (list, np.ndarray)):
        start_positions = [start_positions]
    if not isinstance(end_positions, (list, np.ndarray)):
        end_positions = [end_positions]
    if not isinstance(haplotypes, (list, np.ndarray)):
        haplotypes = [haplotypes]

    start_positions = np.array(start_positions)
    end_positions = np.array(end_positions)
    haplotypes = np.array(haplotypes, dtype=object)

    for haplotype in [1, 2, None]:
        mask = (haplotypes == haplotype)
        if not np.any(mask):
            continue

        starts = start_positions[mask]
        ends = end_positions[mask]

        add_coverage_numba(coverage_arrays[haplotype], starts, ends, coverage_binsize)


def coverage_arrays_to_dict(coverage_arrays, contig):
    result = {contig: {}}
    for haplotype, array in coverage_arrays.items():
        non_zero_indices = np.nonzero(array)[0]
        result[contig][haplotype] = {
            int(idx): int(array[idx]) for idx in non_zero_indices
        }
    return result


def create_coverage_arrays(max_position, coverage_binsize):
    max_bins = (max_position // coverage_binsize) + 1
    return {
        1: np.zeros(max_bins, dtype=np.uint32),
        2: np.zeros(max_bins, dtype=np.uint32),
        None: np.zeros(max_bins, dtype=np.uint32)
    }


def convert_coverage_dict_to_list(coverage_dict, contig_length, coverage_binsize):
    """Convert a coverage dictionary back to a coverage array"""
    coverage_array = {
        1: [0] * ceil(contig_length / coverage_binsize),
        2: [0] * ceil(contig_length / coverage_binsize),
        None: [0] * ceil(contig_length / coverage_binsize)
    }

    max_index = ceil(contig_length / coverage_binsize)

    for haplotype, haplotype_dict in coverage_dict.items():
        if haplotype not in {1, 2, None}:
            print(f"Error: Invalid haplotype {haplotype}. Skipping this haplotype.")
            continue

        for position, coverage in haplotype_dict.items():
            if position < max_index:
                coverage_array[haplotype][position] = coverage
            else:
                print(
                    f"Warning: position {position} exceeds the max_index {max_index} for haplotype {haplotype}. Skipping this position.")

    return coverage_array


def merge_coverage_dicts(coverage_list, contig):
    merged = {1: {}, 2: {}, None: {}}

    for coverage_dict in coverage_list:
        if contig in coverage_dict:
            for haplotype in [1, 2, None]:
                if haplotype in coverage_dict[contig]:
                    for bin_id, count in coverage_dict[contig][haplotype].items():
                        if bin_id in merged[haplotype]:
                            merged[haplotype][bin_id] += count
                        else:
                            merged[haplotype][bin_id] = count

    return merged


@jit(nopython=True)
def accumulate_coverage_numba(bin_ids, counts, target_array):
    for i in range(len(bin_ids)):
        if bin_ids[i] < len(target_array):
            target_array[bin_ids[i]] += counts[i]


def merge_coverage_arrays_efficient(coverage_results, contig, max_position, coverage_binsize):
    max_bins = (max_position // coverage_binsize) + 1
    merged_arrays = {
        1: np.zeros(max_bins, dtype=np.uint32),
        2: np.zeros(max_bins, dtype=np.uint32),
        None: np.zeros(max_bins, dtype=np.uint32)
    }

    for haplotype in [1, 2, None]:
        total_size = sum(
            len(coverage_dict.get(contig, {}).get(haplotype, {}))
            for coverage_dict in coverage_results
            if contig in coverage_dict
        )

        if total_size == 0:
            continue

        all_bin_ids = np.empty(total_size, dtype=np.int32)
        all_counts = np.empty(total_size, dtype=np.uint32)

        idx = 0
        for coverage_dict in coverage_results:
            if contig in coverage_dict and haplotype in coverage_dict[contig]:
                hap_data = coverage_dict[contig][haplotype]
                if hap_data:
                    n_items = len(hap_data)
                    bin_ids = np.fromiter(hap_data.keys(), dtype=np.int32, count=n_items)
                    counts = np.fromiter(hap_data.values(), dtype=np.uint32, count=n_items)

                    all_bin_ids[idx:idx + n_items] = bin_ids
                    all_counts[idx:idx + n_items] = counts
                    idx += n_items

        if idx > 0:
            accumulate_coverage_numba(
                all_bin_ids[:idx],
                all_counts[:idx],
                merged_arrays[haplotype]
            )

    return merged_arrays