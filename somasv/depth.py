from somasv.utils import ensure_int


def compute_depth(breakpoints, shared_cov_arrays, coverage_binsize):
    """
    Use contig coverage to compute local depth and annotate haplotype information.
    """
    for bp in breakpoints:
        bp['phased_local_depths'] = {
            'tumor': [{1: [0, 0], 2: [0, 0], None: [0, 0]} for _ in range(3)],
            'normal': [{1: [0, 0], 2: [0, 0], None: [0, 0]} for _ in range(3)]
        }

        bp['region_total_depths'] = {
            'tumor': [[0, 0] for _ in range(3)],
            'normal': [[0, 0] for _ in range(3)]
        }

        haplotypes = [1, 2, None]

        for sample in bp['phased_local_depths']:
            for i, (chrom, loc) in enumerate([(bp['start_chr'], bp['start_loc']),
                                              (bp['end_chr'], bp['end_loc'])]):
                if chrom not in shared_cov_arrays[sample]:
                    print(f'WARNING: contig {chrom} not in shared_cov_array!')
                    continue

                centre_bin = int(loc // coverage_binsize)

                for hp in haplotypes:
                    if hp not in shared_cov_arrays[sample][chrom]:
                        continue

                    array_length = len(shared_cov_arrays[sample][chrom][hp])
                    centre_bin = min(max(centre_bin, 0), array_length - 1)

                    if centre_bin > 0:
                        depth_before = shared_cov_arrays[sample][chrom][hp][centre_bin - 1]
                        bp['phased_local_depths'][sample][0][hp][i] += depth_before
                        bp['region_total_depths'][sample][0][i] += depth_before

                    depth_current = shared_cov_arrays[sample][chrom][hp][centre_bin]
                    bp['phased_local_depths'][sample][1][hp][i] += depth_current
                    bp['region_total_depths'][sample][1][i] += depth_current

                    if centre_bin < array_length - 1:
                        depth_after = shared_cov_arrays[sample][chrom][hp][centre_bin + 1]
                        bp['phased_local_depths'][sample][2][hp][i] += depth_after
                        bp['region_total_depths'][sample][2][i] += depth_after

        bp['allele_fractions'] = {}
        for sample, regions in bp['phased_local_depths'].items():
            dp_at = regions[1]
            af = [None, None]

            for i in [0, 1]:
                total_depth = bp['region_total_depths'][sample][1][i]
                support_count = bp.get('aln_support_counts', {}).get(sample, 0)

                if total_depth > 0:
                    af[i] = round(support_count / total_depth, 3)
                    af[i] = min(af[i], 1.0)
                else:
                    af[i] = 0.0

            bp['allele_fractions'][sample] = af

        for sample in ['tumor', 'normal']:
            if sample in bp['phased_local_depths']:
                total_hp_at = []
                for haplotype in haplotypes:
                    depth_at_end = bp['phased_local_depths'][sample][1][haplotype][1]
                    total_hp_at.append(str(depth_at_end))

                total_hp_at_str = ",".join(total_hp_at)

                if sample == 'tumor':
                    bp['TUMOR_TOTAL_HP_AT'] = total_hp_at_str
                elif sample == 'normal':
                    bp['NORMAL_TOTAL_HP_AT'] = total_hp_at_str

        for sample in ['tumor', 'normal']:
            if sample in bp['region_total_depths']:
                total_depth_at = [
                    str(bp['region_total_depths'][sample][1][0]),
                    str(bp['region_total_depths'][sample][1][1])
                ]
                total_depth_str = ",".join(total_depth_at)

                if sample == 'tumor':
                    bp['TUMOR_TOTAL_DEPTH_AT'] = total_depth_str
                elif sample == 'normal':
                    bp['NORMAL_TOTAL_DEPTH_AT'] = total_depth_str

        for sample in ['tumor', 'normal']:
            alt_hp_values = ['0', '0', str(bp.get('aln_support_counts', {}).get(sample, 0))]
            alt_hp_str = ",".join(alt_hp_values)

            if sample == 'tumor':
                bp['TUMOR_ALT_HP'] = alt_hp_str
            elif sample == 'normal':
                bp['NORMAL_ALT_HP'] = alt_hp_str

    return breakpoints


def debug_read_counts(final_breakpoints):
    for i, bp in enumerate(final_breakpoints):
        print(f"\n=== Breakpoint {i} ===")
        print(f"total_read_counts: {bp.get('total_read_counts', 'Not found')}")

        start_cluster = bp.get('start_cluster')
        if start_cluster:
            manual_counts = {}
            for cluster_bp in start_cluster['breakpoints']:
                sample = cluster_bp['sample']
                if sample not in manual_counts:
                    manual_counts[sample] = set()
                manual_counts[sample].add(cluster_bp['read_name'])

            manual_totals = {sample: len(reads) for sample, reads in manual_counts.items()}
            print(f"Manual start_cluster counts: {manual_totals}")

        start_stats = bp.get('start_stats_all')
        if start_stats:
            print(f"start_stats_all total_read_counts: {start_stats.get('total_read_counts', 'Not found')}")

        print(
            f"Difference exists: {bp.get('total_read_counts') != start_stats.get('total_read_counts') if start_stats else 'Cannot compare'}")