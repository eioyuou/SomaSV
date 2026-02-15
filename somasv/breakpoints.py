import numpy as np
from math import ceil

from somasv.utils import (
    calculate_statistics_numba, calculate_mean_numba, calculate_median_numba
)


def build_sv_candidate(locations, source, read_name, read_quality, sample, breakpoint_notation, haplotypes=None,
                       phase_sets=None, insert=None):
    """ Create a dictionary representing a breakpoint """
    start_chr = locations[0]['chr']
    start_loc = int(locations[0]['loc'])
    end_chr = locations[1]['chr']
    end_loc = int(locations[1]['loc'])

    is_nearby_event = True if (start_chr == end_chr and abs(start_loc - end_loc) <= 150) else False

    insert_length = len(insert) if insert else 0
    inserted_sequence = insert if insert else None

    haplotypes = haplotypes if haplotypes is not None else [None, None]
    phase_sets = phase_sets if phase_sets is not None else [None, None]

    potential_breakpoint = {
        "start_chr": start_chr,
        "start_loc": start_loc,
        "end_chr": end_chr,
        "end_loc": end_loc,
        "is_nearby_event": is_nearby_event,
        "inserted_sequence": inserted_sequence,
        "insert_length": insert_length,
        "source": source,
        "read_name": read_name,
        "mapq": read_quality,
        "sample": sample,
        "breakpoint_notation": breakpoint_notation,
        "haplotypes": haplotypes if (haplotypes[0] or haplotypes[1]) else None,
        "phase_sets": phase_sets if (phase_sets[0] or phase_sets[1]) else None
    }

    return potential_breakpoint


def built_sv_end_breakpoint(breakpoint):
    """ Reverse breakpoint start and end positions """
    locations = [
        {'chr': breakpoint['start_chr'], 'loc': breakpoint['start_loc']},
        {'chr': breakpoint['end_chr'], 'loc': breakpoint['end_loc']}
    ]

    end_breakpoint = build_sv_candidate(
        locations=locations,
        source=breakpoint['source'],
        read_name=breakpoint['read_name'],
        read_quality=breakpoint['mapq'],
        sample=breakpoint['sample'],
        breakpoint_notation=breakpoint['breakpoint_notation'],
        haplotypes=breakpoint.get('haplotypes', [None, None]),
        phase_sets=breakpoint.get('phase_sets', [None, None]),
        insert=breakpoint.get('inserted_sequence')
    )

    return end_breakpoint


def calculate_comprehensive_cluster_stats(cluster):
    """
    Calculate comprehensive statistics for a cluster
    """
    if not cluster or not cluster.get('breakpoints'):
        return None

    breakpoints = cluster['breakpoints']

    starts = []
    mapqs = []
    event_sizes = []
    samples = set()
    breakpoint_notations = {}
    aln_support_counts = {}
    read_support_counts = {}
    sources = set()
    insert_lengths = []

    unique_reads_per_sample = {}

    haplotype_stats = {1: 0, 2: 0, None: 0}
    phase_sets = set()

    qualities = []

    for bp in breakpoints:
        starts.append(bp['start_loc'])
        mapqs.append(bp.get('mapq', 0))
        qualities.append(bp.get('mapq', 0))

        sample = bp.get('sample', 'unknown')
        samples.add(sample)

        notation = bp.get('breakpoint_notation', 'unknown')
        breakpoint_notations[notation] = breakpoint_notations.get(notation, 0) + 1

        read_name = bp.get('read_name')
        if read_name:
            if sample not in read_support_counts:
                read_support_counts[sample] = set()
            read_support_counts[sample].add(read_name)

            if sample not in unique_reads_per_sample:
                unique_reads_per_sample[sample] = set()
            unique_reads_per_sample[sample].add(read_name)

        aln_support_counts[sample] = aln_support_counts.get(sample, 0) + 1

        source = bp.get('source', 'unknown')
        sources.add(source)

        if notation == "<INS>":
            insert_len = bp.get('insert_length', 0)
            event_sizes.append(insert_len)
            insert_lengths.append(insert_len)
        elif notation in ["+", "-"]:
            insert_len = bp.get('insert_length', 0)
            event_sizes.append(insert_len)
        elif bp.get('start_chr') == bp.get('end_chr'):
            if bp.get('end_loc') is not None:
                sv_size = abs(bp['start_loc'] - bp['end_loc'])
                event_sizes.append(sv_size)
        else:
            event_sizes.append(0)

        haplotypes = bp.get('haplotypes', [None, None])
        if haplotypes and len(haplotypes) >= 2:
            for hp in haplotypes:
                if hp in haplotype_stats:
                    haplotype_stats[hp] += 1

        phase_sets_bp = bp.get('phase_sets', [None, None])
        if phase_sets_bp:
            for ps in phase_sets_bp:
                if ps is not None:
                    phase_sets.add(ps)

    total_read_counts = {}
    for sample, reads in unique_reads_per_sample.items():
        total_read_counts[sample] = len(reads)

    read_counts_by_sample = {sample: len(reads) for sample, reads in read_support_counts.items()}

    if starts and mapqs and event_sizes:
        starts_std, mapq_mean, event_size_std, event_size_median, event_size_mean = \
            calculate_statistics_numba(starts, mapqs, event_sizes)
    else:
        starts_std = mapq_mean = event_size_std = event_size_median = event_size_mean = 0.0

    stats = {
        'total_breakpoints': len(breakpoints),
        'unique_samples': len(samples),
        'sample_list': sorted(list(samples)),
        'unique_sources': len(sources),
        'source_list': sorted(list(sources)),
        'starts_std_dev': round(starts_std, 3),
        'starts_mean': round(calculate_mean_numba(starts), 3) if starts else 0.0,
        'starts_median': round(calculate_median_numba(starts), 3) if starts else 0.0,
        'starts_range': [min(starts), max(starts)] if starts else [0, 0],
        'mapq_mean': round(mapq_mean, 3),
        'mapq_median': round(calculate_median_numba(mapqs), 3) if mapqs else 0.0,
        'mapq_std_dev': round(np.std(mapqs), 3) if mapqs else 0.0,
        'mapq_range': [min(mapqs), max(mapqs)] if mapqs else [0, 0],
        'event_size_mean': round(event_size_mean, 3),
        'event_size_median': round(event_size_median, 3),
        'event_size_std_dev': round(event_size_std, 3),
        'event_size_range': [min(event_sizes), max(event_sizes)] if event_sizes else [0, 0],
        'insert_length_stats': {
            'count': len(insert_lengths),
            'mean': round(calculate_mean_numba(insert_lengths), 3) if insert_lengths else 0.0,
            'median': round(calculate_median_numba(insert_lengths), 3) if insert_lengths else 0.0,
            'std_dev': round(np.std(insert_lengths), 3) if insert_lengths else 0.0,
            'range': [min(insert_lengths), max(insert_lengths)] if insert_lengths else [0, 0]
        },
        'breakpoint_notation_counts': breakpoint_notations,
        'dominant_notation': max(breakpoint_notations.items(), key=lambda x: x[1]) if breakpoint_notations else None,
        'aln_support_counts': aln_support_counts,
        'read_support_counts': read_counts_by_sample,
        'total_read_counts': total_read_counts,
        'total_alignments': sum(aln_support_counts.values()),
        'total_unique_reads': sum(total_read_counts.values()),
        'haplotype_distribution': haplotype_stats,
        'phase_sets_count': len(phase_sets),
        'phase_sets_list': sorted(list(phase_sets)) if phase_sets else [],
        'cluster_span': cluster.get('end', 0) - cluster.get('start', 0),
        'cluster_chr': cluster.get('chr', 'unknown'),
        'cluster_start': cluster.get('start', 0),
        'cluster_end': cluster.get('end', 0),
        'breakpoint_density': len(breakpoints) / max(1, cluster.get('end', 1) - cluster.get('start', 0)) * 1000,
        'max_sample_support': max(read_counts_by_sample.values()) if read_counts_by_sample else 0,
        'min_sample_support': min(read_counts_by_sample.values()) if read_counts_by_sample else 0,
        'support_ratio': max(read_counts_by_sample.values()) / len(breakpoints) if read_counts_by_sample and breakpoints else 0.0
    }

    return stats


def update_consensus_breakpoint_with_cluster_stats(consensus_breakpoint):
    """
    Update the consensus breakpoint by adding start_stats_all and end_stats_all
    """
    start_cluster = consensus_breakpoint.get('start_cluster')
    end_cluster = consensus_breakpoint.get('end_cluster')

    if start_cluster:
        start_stats = calculate_comprehensive_cluster_stats(start_cluster)
        consensus_breakpoint['start_stats_all'] = start_stats
    else:
        consensus_breakpoint['start_stats_all'] = None

    if end_cluster and end_cluster != start_cluster:
        end_stats = calculate_comprehensive_cluster_stats(end_cluster)
        consensus_breakpoint['end_stats_all'] = end_stats
    else:
        consensus_breakpoint['end_stats_all'] = consensus_breakpoint['start_stats_all']

    return consensus_breakpoint


def create_consensus_breakpoint_with_stats(locations, source, start_cluster, end_cluster, counts, breakpoint_notation,
                                           inserts=None):
    """
    Create a consensus breakpoint with full statistical information.
    """
    from somasv.clustering import calculate_cluster_statistics

    start_chr = locations[0]['chr']
    start_loc = int(locations[0]['loc'])
    end_chr = locations[1]['chr']
    end_loc = int(locations[1]['loc'])

    calculate_cluster_statistics(start_cluster)
    if end_cluster and end_cluster != start_cluster:
        calculate_cluster_statistics(end_cluster)

    read_support_counts = {'normal': 0, 'tumor': 0}
    supporting_reads = counts[0]
    for sample, reads in supporting_reads.items():
        read_support_counts[sample] += len(reads)

    sv_length = None
    if breakpoint_notation == "<INS>":
        sv_length = str(int(calculate_mean_numba([
            bp['insert_length'] for bp in start_cluster['breakpoints']
            if bp['breakpoint_notation'] == "<INS>"
        ])))
    elif breakpoint_notation in ["-", "+"]:
        sv_length = 0
    elif start_chr != end_chr:
        sv_length = 0
    else:
        sv_length = abs(int(end_loc) - int(start_loc))

    consensus_breakpoint = {
        "start_chr": start_chr,
        "start_loc": start_loc,
        "end_chr": end_chr,
        "end_loc": end_loc,
        "source": source,
        "breakpoint_notation": breakpoint_notation,
        "inserted_sequences": inserts,
        "start_cluster": start_cluster,
        "end_cluster": end_cluster if end_cluster else start_cluster,
        "supporting_reads": supporting_reads,
        "aln_support_counts": counts[1],
        "phasing": counts[2],
        "total_read_counts": counts[3],
        "read_support_counts": read_support_counts,
        "sv_length": sv_length,
        "count": None,
        "allele_fractions": {}
    }

    consensus_breakpoint = update_consensus_breakpoint_with_cluster_stats(consensus_breakpoint)

    return consensus_breakpoint


def create_consensus_breakpoint(locations, source, start_cluster, end_cluster, counts, breakpoint_notation,
                                inserts=None):
    return create_consensus_breakpoint_with_stats(locations, source, start_cluster, end_cluster, counts,
                                                  breakpoint_notation, inserts)


def batch_update_consensus_breakpoints_with_stats(consensus_breakpoints):
    updated_breakpoints = []
    for bp in consensus_breakpoints:
        updated_bp = update_consensus_breakpoint_with_cluster_stats(bp)
        updated_breakpoints.append(updated_bp)
    return updated_breakpoints


def format_cluster_stats_summary(stats):
    if not stats:
        return "No statistics available"

    summary = f"""
    Cluster Statistics Summary:
    ========================
    Total Breakpoints: {stats['total_breakpoints']}
    Samples: {', '.join(stats['sample_list'])} (n={stats['unique_samples']})
    Sources: {', '.join(stats['source_list'])} (n={stats['unique_sources']})

    Position Stats:
    - Range: {stats['starts_range'][0]:,} - {stats['starts_range'][1]:,} bp
    - Mean ± SD: {stats['starts_mean']:,.1f} ± {stats['starts_std_dev']:.1f} bp
    - Span: {stats['cluster_span']:,} bp

    Quality Stats:
    - MAPQ Mean ± SD: {stats['mapq_mean']:.1f} ± {stats['mapq_std_dev']:.1f}
    - MAPQ Range: {stats['mapq_range'][0]} - {stats['mapq_range'][1]}

    Event Size Stats:
    - Mean ± SD: {stats['event_size_mean']:,.1f} ± {stats['event_size_std_dev']:.1f} bp
    - Median: {stats['event_size_median']:,.1f} bp

    Breakpoint Types:
    """

    for notation, count in stats['breakpoint_notation_counts'].items():
        summary += f"    {notation}: {count}\n"

    summary += f"""
    Support:
    - Total Alignments: {stats['total_alignments']}
    - Unique Reads: {stats['total_unique_reads']}
    - Density: {stats['breakpoint_density']:.2f} breakpoints/kb
    """

    return summary