import multiprocessing

from somasv.utils import calculate_median_numba
from somasv.clustering import (
    initialize_cluster, merge_breakpoint_into_cluster, cluster_by_insert_length,
    group_related_end_breakpoints
)
from somasv.breakpoints import built_sv_end_breakpoint, create_consensus_breakpoint


def count_samples_and_phase(source_breakpoints):
    """ Given a set of breakpoints, return count of labels, all alignments, and their phasing """
    aln_counts = {}
    read_counts = {}
    seen_reads = set()
    start_phase, end_phase = {}, {}

    for bp in source_breakpoints:
        if bp['sample'] not in start_phase:
            start_phase[bp['sample']] = {'HP': {1: 0, 2: 0, None: 0}, 'PS': set()}
        if bp['sample'] not in end_phase:
            end_phase[bp['sample']] = {'HP': {1: 0, 2: 0, None: 0}, 'PS': set()}

        if bp['haplotypes']:
            start_phase[bp['sample']]['HP'][bp['haplotypes'][0]] += 1
            if bp['phase_sets'][0] is not None:
                start_phase[bp['sample']]['PS'].add(bp['phase_sets'][0])

            end_phase[bp['sample']]['HP'][bp['haplotypes'][1]] += 1
            if bp['phase_sets'][1] is not None:
                end_phase[bp['sample']]['PS'].add(bp['phase_sets'][1])
        else:
            start_phase[bp['sample']]['HP'][None] += 1
            end_phase[bp['sample']]['HP'][None] += 1

        aln_counts[bp['sample']] = aln_counts.get(bp['sample'], 0) + 1

        read_key = (bp['read_name'], bp['sample'])
        if read_key not in seen_reads:
            read_counts.setdefault(bp['sample'], []).append(bp['read_name'])
            seen_reads.add(read_key)

    return read_counts, aln_counts, [start_phase, end_phase]


def process_insertion_like_breakpoints(cluster, insertion_like_breakpoints, total_read_counts, min_support):
    """ Process insertion/single-end breakpoints """
    breakpoints_for_end_chrom = []

    if len(insertion_like_breakpoints) >= min_support:
        insertion_like_clusters = cluster_by_insert_length(insertion_like_breakpoints, 0.75)

        for ins_cluster in insertion_like_clusters:
            ins_read_counts, ins_aln_counts, ins_phasing_counts = count_samples_and_phase(ins_cluster)

            if max(len(count) for count in ins_read_counts.values()) >= min_support:
                num_insertions = sum(b['breakpoint_notation'] == "<INS>" for b in ins_cluster)

                if num_insertions >= 2:
                    event_info = {'starts': [], 'inserts': [], 'sources': {}}
                    start_cluster = None

                    for bp in ins_cluster:
                        event_info['starts'].append(bp['start_loc'])
                        event_info['inserts'].append(bp['inserted_sequence'])
                        event_info['sources'].setdefault(bp['source'], True)

                        if not start_cluster:
                            start_cluster = initialize_cluster(bp)
                        else:
                            merge_breakpoint_into_cluster(start_cluster, bp)

                    median_start = calculate_median_numba(event_info['starts'])
                    consensus_source = "/".join(sorted(event_info['sources'].keys()))

                    breakpoints_for_end_chrom.append(create_consensus_breakpoint(
                        [{'chr': cluster['chr'], 'loc': median_start}, {'chr': cluster['chr'], 'loc': median_start}],
                        consensus_source, start_cluster, None,
                        [ins_read_counts, ins_aln_counts, ins_phasing_counts, total_read_counts],
                        "<INS>", event_info['inserts']
                    ))
                else:
                    sorted_breakpoints = {}
                    for bp in ins_cluster:
                        sorted_breakpoints.setdefault(bp['breakpoint_notation'], []).append(bp)

                    for notation_type, breakpoints in sorted_breakpoints.items():
                        event_info = {'starts': [], 'inserts': [], 'sources': {}}
                        start_cluster = None

                        for bp in breakpoints:
                            event_info['starts'].append(bp['start_loc'])
                            event_info['inserts'].append(bp['inserted_sequence'])
                            event_info['sources'].setdefault(bp['source'], True)

                            if not start_cluster:
                                start_cluster = initialize_cluster(bp)
                            else:
                                merge_breakpoint_into_cluster(start_cluster, bp)

                        median_start = calculate_median_numba(event_info['starts'])
                        consensus_source = "/".join(sorted(event_info['sources'].keys()))
                        sbnd_read_counts, sbnd_aln_counts, sbnd_phasing_counts = count_samples_and_phase(breakpoints)

                        if max(len(count) for count in sbnd_read_counts.values()) >= min_support:
                            breakpoints_for_end_chrom.append(create_consensus_breakpoint(
                                [{'chr': cluster['chr'], 'loc': median_start},
                                 {'chr': cluster['chr'], 'loc': median_start}],
                                consensus_source, start_cluster, None,
                                [sbnd_read_counts, sbnd_aln_counts, sbnd_phasing_counts, total_read_counts],
                                notation_type, event_info['inserts']
                            ))

    return breakpoints_for_end_chrom


def process_end_breakpoints(end_chrom, end_chrom_breakpoints, end_extension, total_read_counts, tumor_bamfile,
                            normal_bamfile, min_support, min_length):
    """ Process reverse breakpoints """
    breakpoints_for_end_chrom = []

    end_breakpoints = [built_sv_end_breakpoint(b) for b in end_chrom_breakpoints if
                       b['breakpoint_notation'] not in ["<INS>", "+", "-"]]
    end_breakpoints.sort(key=lambda x: x['end_loc'])

    _, cluster_stack = group_related_end_breakpoints(end_chrom, end_breakpoints, tumor_bamfile, normal_bamfile,
                                                     end_extension)

    for end_chrom_cluster in cluster_stack:
        sorted_breakpoints = {}
        for bp in end_chrom_cluster['breakpoints']:
            sorted_breakpoints.setdefault(bp['breakpoint_notation'], []).append(bp)

        for bp_type, breakpoints in sorted_breakpoints.items():
            read_counts, aln_counts, phasing_counts = count_samples_and_phase(breakpoints)

            if max([len(count) for count in read_counts.values()]) >= min_support:
                start_cluster = None
                source = {}

                for bp in breakpoints:
                    if not start_cluster:
                        start_cluster = initialize_cluster(bp)
                    else:
                        merge_breakpoint_into_cluster(start_cluster, bp)
                    source.setdefault(bp['source'], True)

                median_start = calculate_median_numba([bp['start_loc'] for bp in breakpoints])
                median_end = calculate_median_numba([bp['end_loc'] for bp in breakpoints])
                consensus_source = "/".join(sorted(source.keys()))

                new_breakpoint = create_consensus_breakpoint(
                    [{'chr': start_cluster['chr'], 'loc': median_start},
                     {'chr': end_chrom_cluster['chr'], 'loc': median_end}],
                    consensus_source, start_cluster, end_chrom_cluster,
                    [read_counts, aln_counts, phasing_counts, total_read_counts],
                    bp_type
                )

                if new_breakpoint['sv_length'] >= min_length or (
                        new_breakpoint['start_chr'] != new_breakpoint['end_chr'] and new_breakpoint['sv_length'] == 0):
                    breakpoints_for_end_chrom.append(new_breakpoint)

    return breakpoints_for_end_chrom


def svtype_breakpoints(breakpoints_for_end_chrom):
    """ Genotype determination """
    final_breakpoints = []
    counts = {}

    for new_breakpoint in breakpoints_for_end_chrom:
        counts[new_breakpoint['breakpoint_notation']] = counts.get(new_breakpoint['breakpoint_notation'], 0) + 1

    bp_types = counts.keys()
    if "-" in bp_types or "+" in bp_types:
        major_types = [t for t in bp_types if t not in ("+", "-")]
        if major_types:
            major_types_bps = [b for b in breakpoints_for_end_chrom if b['breakpoint_notation'] not in ("+", "-")]
            final_breakpoints.extend(major_types_bps)
        else:
            final_breakpoints.extend(breakpoints_for_end_chrom)
    else:
        final_breakpoints.extend(breakpoints_for_end_chrom)

    return final_breakpoints


def call_breakpoints(clusters, end_extension, tumor_bamfile, normal_bamfile, min_length, min_support, chrom):
    """ Identify consensus breakpoints from clusters """
    final_breakpoints = []

    for cluster in clusters:
        breakpoints_by_end_chrom = {}
        total_read_counts = {}
        unique_reads_per_sample = {}

        for bp in cluster['breakpoints']:
            breakpoints_by_end_chrom.setdefault(bp['end_chr'], []).append(bp)

            if bp['sample'] not in unique_reads_per_sample:
                unique_reads_per_sample[bp['sample']] = set()
            unique_reads_per_sample[bp['sample']].add(bp['read_name'])

        for sample, reads in unique_reads_per_sample.items():
            total_read_counts[sample] = len(reads)

        for end_chrom, end_chrom_breakpoints in breakpoints_by_end_chrom.items():
            breakpoints_for_end_chrom = []

            if len(end_chrom_breakpoints) >= min_support:
                insertion_like_breakpoints = [b for b in end_chrom_breakpoints if
                                              b['breakpoint_notation'] in ["<INS>", "+", "-"]]
                breakpoints_for_end_chrom.extend(
                    process_insertion_like_breakpoints(cluster, insertion_like_breakpoints, total_read_counts,
                                                       min_support))

                breakpoints_for_end_chrom.extend(
                    process_end_breakpoints(end_chrom, end_chrom_breakpoints, end_extension, total_read_counts,
                                            tumor_bamfile, normal_bamfile,
                                            min_support, min_length))

            final_breakpoints.extend(svtype_breakpoints(breakpoints_for_end_chrom))

    breakpoints_call_set = {}
    breakpoints_call_set.setdefault(chrom, []).extend(final_breakpoints)
    print(f"Chromosome {chrom} identified {len(final_breakpoints)} breakpoints.")
    return breakpoints_call_set


def call_breakpoints_task(args):
    """Multiprocess task function for calling call_breakpoints"""
    print("Processing chromosome:", args[4])
    clusters, end_extension, min_length, min_support, chrom, tumor_bamfile, normal_bamfile = args
    return call_breakpoints(clusters, end_extension, tumor_bamfile, normal_bamfile, min_length, min_support, chrom)


def parallel_call_breakpoints(clusters_by_chrom, tumor_bamfile, normal_bamfile, end_extension, min_length, min_support,
                              num_processes=None):
    tasks = [(clusters, end_extension, min_length, min_support, chrom, tumor_bamfile, normal_bamfile) for
             chrom, clusters in clusters_by_chrom.items()]

    with multiprocessing.Pool(processes=num_processes) as pool:
        results = pool.map(call_breakpoints_task, tasks)

    final_breakpoints_by_chrom = {}
    for result in results:
        for chrom, breakpoints in result.items():
            if chrom not in final_breakpoints_by_chrom:
                final_breakpoints_by_chrom[chrom] = breakpoints
            else:
                final_breakpoints_by_chrom[chrom].extend(breakpoints)

    return final_breakpoints_by_chrom