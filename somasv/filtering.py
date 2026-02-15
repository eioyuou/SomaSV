import math
import numpy as np

from somasv.utils import (
    ensure_float, ensure_int,
    get_breakpoint_af, get_breakpoint_support, get_breakpoint_mapq,
    get_breakpoint_clustered_reads, get_breakpoint_std_dev
)
from somasv.vcf_io import (
    load_bed_regions, get_support_read,
    generate_vcf_header, generate_vcf_line, write_somatic_vcf
)


def filter_variants(variant, coverage):
    svlen = abs(int(variant.get('sv_length', 0)))
    if svlen < 50 and svlen > 0:
        return False

    tumor_af = variant.get('allele_fractions', {}).get('tumor', [0, 0])
    normal_af = variant.get('allele_fractions', {}).get('normal', [0, 0])

    tumor_af_max = max(tumor_af[0], tumor_af[1]) if isinstance(tumor_af, list) else tumor_af
    normal_af_max = max(normal_af[0], normal_af[1]) if isinstance(normal_af, list) else normal_af

    min_support = math.ceil(coverage * 0.1) + 1
    tumor_support = variant['read_support_counts'].get('tumor', 0)
    normal_support = variant['read_support_counts'].get('normal', 0)

    if normal_af_max > 0.1:
        return False
    if tumor_af_max < 0.1:
        return False
    if variant['start_cluster']['stats']['mapq_mean'] < 40.0:
        return False
    if 'end_cluster' in variant and variant['end_cluster']['stats'].get('mapq_mean', 0) < 40.0:
        return False
    if normal_support > 1:
        return False
    if tumor_support < min_support:
        return False
    if variant['total_read_counts'].get('normal', 0) > 1:
        return False

    return True


def filter_somatic_af(bps):
    afs = []
    for bp in bps:
        if 'allele_fractions' not in bp:
            continue

        tumor_afs = bp['allele_fractions'].get('tumor', [0, 0])
        normal_afs = bp['allele_fractions'].get('normal', [0, 0])
        tumor_af = np.mean(tumor_afs)
        normal_af = np.mean(normal_afs)

        if tumor_af >= 0.1 and normal_af <= 0.01:
            afs.append(tumor_af)

    return np.mean(afs) if afs else 0


def filter_sv_length(bp, min_length=50):
    svlen = abs(ensure_int(bp.get('sv_length', 0)))
    if svlen < min_length and svlen > 0:
        return False
    return True


def depth_all(bps):
    depth = []
    for bp in bps:
        depth.append(min(bp['region_total_depths']['tumor'][1]))
    return np.mean(depth) * 0.1


def save_stats_to_txt(global_stats, mean_af, output_file):
    with open(output_file, 'w') as f:
        f.write("=== Global Statistics ===\n")
        for key, value in global_stats.items():
            f.write(f"{key}: {value}\n")
        f.write("\n=== Additional Metrics ===\n")
        f.write(f"mean_af: {mean_af:.5f}\n")


class SomaticFilter:
    def __init__(self, platform, mean_af, tumor_global_coverage, normal_global_coverage, tumor_mean_af, normal_mean_af):
        self.platform = platform
        self.mean_af = mean_af
        self.tumor_global_coverage = tumor_global_coverage
        self.normal_global_coverage = normal_global_coverage
        self.tumor_mean_af = tumor_mean_af
        self.normal_mean_af = normal_mean_af
        self.use_strict_mode = self._should_use_strict_mode()

    def _should_use_strict_mode(self):
        if self.mean_af > 0.25:
            return True
        else:
            return False

    def filter_breakpoint(self, bp, bed_support_info=None):
        if not self._basic_prefilter(bp):
            return False, "Failed basic prefilter"
        if bed_support_info:
            return self._filter_bed_region_lenient(bp, bed_support_info)
        else:
            if self.use_strict_mode:
                return self._filter_strict_region_high(bp)
            else:
                return self._filter_strict_region_low(bp)

    def _basic_prefilter(self, bp):
        if not filter_sv_length(bp):
            return False
        if 'allele_fractions' not in bp:
            return False
        return True

    def _filter_strict_region_high(self, bp):
        tumor_af, normal_af = get_breakpoint_af(bp)
        tumor_support, normal_support, support_ratio = get_breakpoint_support(bp)
        origin_mapq, end_mapq = get_breakpoint_mapq(bp)
        clustered_tumor, clustered_normal = get_breakpoint_clustered_reads(bp)

        tumor_coverage = min(ensure_int(bp['region_total_depths']['tumor'][1][0]),
                             ensure_int(bp['region_total_depths']['tumor'][1][1]))
        normal_coverage = min(ensure_int(bp['region_total_depths']['normal'][1][0]),
                              ensure_int(bp['region_total_depths']['normal'][1][1]))

        normal_coverage_ratio = normal_coverage / self.normal_global_coverage
        if normal_coverage_ratio < 0.3 or normal_coverage_ratio > 3.0:
            return False, "abnormal_coverage_ratio"
        if normal_coverage < 2:
            return False, "low_normal_coverage"

        min_support_tumor = math.ceil(tumor_coverage * 0.1) + 1 if tumor_coverage else 3
        max_support_normal = math.floor(normal_coverage * 0.1) if normal_coverage else 1

        if normal_af > 0.05:
            return False, "high_normal_af"
        if tumor_af < 0.05:
            return False, "low_tumor_af"
        if origin_mapq < 55.0 or end_mapq < 55.0:
            return False, "low_mapq"
        if abs(origin_mapq - end_mapq) > 5:
            return False, "mapq_diff_too_large"
        if normal_support > max_support_normal:
            return False, "high_normal_support"
        if tumor_support < min_support_tumor:
            return False, "low_tumor_support"

        try:
            origin_starts_std_dev = bp['start_cluster']['stats']['starts_std_dev'] / (
                    bp['end_cluster']['stats']['event_size_mean'] + 1)
            end_starts_std_dev = bp['end_cluster']['stats']['starts_std_dev'] / (
                    bp['end_cluster']['stats']['event_size_mean'] + 1)
            if origin_starts_std_dev > 10 or end_starts_std_dev > 40:
                return False, "high_std_dev"
        except:
            pass

        if clustered_normal > 1:
            if clustered_normal / (clustered_tumor + 1e-8) > 0.1:
                return False, "high_clustered_normal_ratio"

        min_support_ratio_threshold = 0.15
        min_support_ratio = self.tumor_global_coverage * min_support_ratio_threshold
        if support_ratio < min_support_ratio:
            return False, "low_support_ratio"

        if bp['breakpoint_notation'] == '<INS>':
            starts_std_dev = bp.get('start_stats_all', {}).get('starts_std_dev', 0)
            if starts_std_dev > 30:
                return False, "insertion_high_std_dev"

        return True, "passed_strict_high"

    def _filter_strict_region_low(self, bp):
        tumor_af, normal_af = get_breakpoint_af(bp)
        tumor_support, normal_support, support_ratio = get_breakpoint_support(bp)
        origin_mapq, end_mapq = get_breakpoint_mapq(bp)
        clustered_tumor, clustered_normal = get_breakpoint_clustered_reads(bp)

        normal_coverage = min(ensure_int(bp['region_total_depths']['normal'][1][0]),
                              ensure_int(bp['region_total_depths']['normal'][1][1]))

        if tumor_af < 0.05:
            return False, "low_tumor_af"
        if normal_af > 0.05:
            return False, "high_normal_af"
        if origin_mapq < 55.0 or end_mapq < 55.0:
            return False, "low_mapq"
        if normal_coverage < 2:
            return False, "low_normal_coverage"
        if abs(origin_mapq - end_mapq) > 5:
            return False, "mapq_diff_too_large"

        try:
            origin_starts_std_dev = bp['start_cluster']['stats']['starts_std_dev'] / (
                    bp['end_cluster']['stats']['event_size_mean'] + 1)
            end_starts_std_dev = bp['end_cluster']['stats']['starts_std_dev'] / (
                    bp['end_cluster']['stats']['event_size_mean'] + 1)
            if origin_starts_std_dev > 10 or end_starts_std_dev > 40:
                return False, "high_std_dev"
        except:
            pass

        if clustered_normal > 1:
            if clustered_normal / (clustered_tumor + 1e-8) > 0.1:
                return False, "high_clustered_normal_ratio"

        return True, "passed_strict_low"

    def _filter_bed_region_lenient(self, bp, bed_support_info):
        tumor_af, normal_af = get_breakpoint_af(bp)
        tumor_support, normal_support, support_ratio = get_breakpoint_support(bp)
        origin_mapq, end_mapq = get_breakpoint_mapq(bp)
        clustered_tumor, clustered_normal = get_breakpoint_clustered_reads(bp)

        tumor_coverage = min(ensure_int(bp['region_total_depths']['tumor'][1][0]),
                             ensure_int(bp['region_total_depths']['tumor'][1][1]))

        min_support = math.ceil(tumor_coverage * 0.1) + 1 if tumor_coverage else 3

        if normal_af > 0.05:
            if (normal_af <= 0.07 and tumor_af >= 0.2 and
                    tumor_support >= min_support * 3 and support_ratio >= 10):
                pass
            else:
                return False, "bed_high_normal_af"

        if tumor_af < 0.1:
            if (tumor_af >= 0.08 and tumor_support >= min_support * 3 and normal_support == 0):
                pass
            else:
                return False, "bed_low_tumor_af"

        if origin_mapq < 40.0 or end_mapq < 40.0:
            min_mapq = min(origin_mapq, end_mapq)
            if (min_mapq >= 35.0 and tumor_support >= min_support * 3 and
                    tumor_af >= 0.2 and normal_support == 0):
                pass
            else:
                return False, "bed_low_mapq"

        if normal_support > 0:
            if (normal_support <= 1 and tumor_support >= min_support * 3 and
                    tumor_af >= 0.2 and support_ratio >= min_support):
                pass
            else:
                return False, "bed_high_normal_support"

        if tumor_support < min_support:
            if (tumor_support >= min_support * 0.8 and tumor_af >= 0.3 and normal_support == 0):
                pass
            else:
                return False, "bed_low_tumor_support"

        if clustered_normal > 3:
            if (clustered_normal <= 4 and tumor_support >= min_support * 3 and
                    tumor_af >= 0.3 and normal_support == 0):
                pass
            else:
                return False, "bed_high_clustered_normal"

        return True, "passed_bed_lenient"


def calculate_interpretable_features(bp, tumor_coverage, normal_coverage, platform):
    features = {}

    tumor_af, normal_af = get_breakpoint_af(bp)
    tumor_support, normal_support, support_ratio = get_breakpoint_support(bp)
    origin_mapq, end_mapq = get_breakpoint_mapq(bp)

    features['tumor_af'] = tumor_af
    features['normal_af'] = normal_af
    features['tumor_support'] = tumor_support
    features['normal_support'] = normal_support
    features['support_ratio'] = support_ratio
    features['min_mapq'] = min(origin_mapq, end_mapq)
    features['mapq_diff'] = abs(origin_mapq - end_mapq)

    features['tumor_coverage'] = tumor_coverage
    features['normal_coverage'] = normal_coverage
    features['support_density'] = tumor_support / (tumor_coverage + 1)

    expected_support = max(2, tumor_coverage * 0.06)
    features['support_ratio_expected'] = tumor_support / expected_support
    features['af_support_consistency'] = tumor_af * tumor_coverage / (tumor_support + 1)

    clustered_tumor, clustered_normal = get_breakpoint_clustered_reads(bp)
    features['clustered_normal'] = clustered_normal
    features['normal_noise_ratio'] = clustered_normal / (clustered_tumor + 1)

    if platform == 'ONT':
        try:
            origin_starts_std_dev = bp.get('start_cluster', {}).get('stats', {}).get('starts_std_dev', 0)
            end_starts_std_dev = bp.get('end_cluster', {}).get('stats', {}).get('starts_std_dev', 0)
            event_size_mean = bp.get('end_cluster', {}).get('stats', {}).get('event_size_mean', 1)

            features['std_dev_ratio'] = max(
                origin_starts_std_dev / (event_size_mean + 1),
                end_starts_std_dev / (event_size_mean + 1)
            )
        except:
            features['std_dev_ratio'] = 0

    return features


def calculate_relative_confidence(features, reference_stats=None):
    confidence_factors = {}
    confidence_factors['af_confidence'] = min(1.0, features['tumor_af'] / 0.05)
    confidence_factors['support_confidence'] = min(1.0, features['support_ratio_expected'])
    confidence_factors['quality_confidence'] = min(1.0, features['min_mapq'] / 60.0)
    confidence_factors['noise_confidence'] = 1.0 / (1.0 + features['normal_noise_ratio'] * 10)
    normal_af_penalty = max(0, features['normal_af'] - 0.01) * 20
    confidence_factors['normal_af_confidence'] = max(0.1, 1.0 - normal_af_penalty)

    overall_confidence = sum(confidence_factors.values()) / (len(confidence_factors) + 1)

    return overall_confidence, confidence_factors


def calculate_reference_stats(breakpoints):
    afs = []
    supports = []
    mapqs = []

    for bp in breakpoints:
        if 'allele_fractions' in bp:
            tumor_af, _ = get_breakpoint_af(bp)
            tumor_support, _, _ = get_breakpoint_support(bp)
            origin_mapq, end_mapq = get_breakpoint_mapq(bp)

            afs.append(tumor_af)
            supports.append(tumor_support)
            mapqs.append(min(origin_mapq, end_mapq))

    return {
        'af_percentiles': np.percentile(afs, [10, 25, 50, 75, 90]) if afs else [0] * 5,
        'support_median': np.median(supports) if supports else 0,
        'mapq_median': np.median(mapqs) if mapqs else 0,
        'total_variants': len(breakpoints)
    }


def should_use_strict_mode(reference_stats, global_stats):
    if not global_stats:
        return False
    high_coverage = global_stats.get('tumor_global_coverage', 0) > 30
    sufficient_variants = reference_stats['total_variants'] > 100
    good_quality = reference_stats['mapq_median'] > 50
    return high_coverage and sufficient_variants and good_quality


def compare_af_distributions(bps):
    tumor_afs = []
    normal_afs = []

    for bp in bps:
        if 'allele_fractions' in bp:
            tumor_af_list = bp['allele_fractions'].get('tumor', [])
            normal_af_list = bp['allele_fractions'].get('normal', [])

            if tumor_af_list:
                tumor_af = max(tumor_af_list)
                if 0.1 <= tumor_af <= 0.6:
                    tumor_afs.append(tumor_af)

            if normal_af_list:
                normal_af = max(normal_af_list)
                if 0.1 <= normal_af <= 0.6:
                    normal_afs.append(normal_af)

    tumor_mean_af = sum(tumor_afs) / len(tumor_afs) if tumor_afs else None
    normal_mean_af = sum(normal_afs) / len(normal_afs) if normal_afs else None

    return tumor_mean_af, normal_mean_af


def output_filtering_stats(passed_breakpoints, reference_stats, outdir):
    stats_file = f"{outdir}.filtering_stats_data.txt"

    with open(stats_file, 'w') as f:
        f.write("=== Somatic Filtering Statistics ===\n")
        f.write(f"Reference AF percentiles: {reference_stats['af_percentiles']}\n")
        f.write(f"Reference support median: {reference_stats['support_median']}\n")
        f.write(f"Reference MAPQ median: {reference_stats['mapq_median']}\n")
        f.write(f"Total variants processed: {reference_stats['total_variants']}\n")
        f.write(f"Variants passed: {len(passed_breakpoints)}\n")

        if passed_breakpoints:
            f.write("\n=== Confidence Factor Distribution ===\n")
            for factor_name in ['af_confidence', 'support_confidence', 'quality_confidence', 'noise_confidence']:
                values = [bp['CONFIDENCE_FACTORS'][factor_name] for bp in passed_breakpoints if
                          'CONFIDENCE_FACTORS' in bp]
                if values:
                    f.write(f"{factor_name}: mean={np.mean(values):.3f}, median={np.median(values):.3f}\n")


def filter_somatic_variants_improved(breakpoints, outdir, platform, ref_fasta_path,
                                     global_stats=None, bed_file=None):
    reference_stats = calculate_reference_stats(breakpoints)

    mean_af = filter_somatic_af(breakpoints)
    tumor_mean_af, normal_mean_af = compare_af_distributions(breakpoints)
    print(f"Mean AF for somatic variants: {mean_af:.3f}")
    tumor_global_coverage = global_stats.get('tumor_global_coverage', 30)
    normal_global_coverage = global_stats.get('normal_global_coverage', 30)
    save_stats_to_txt(global_stats, mean_af, outdir + ".global_stats.txt")

    somatic_filter = SomaticFilter(
        platform=platform,
        mean_af=mean_af,
        tumor_global_coverage=tumor_global_coverage,
        normal_global_coverage=normal_global_coverage,
        tumor_mean_af=tumor_mean_af,
        normal_mean_af=normal_mean_af,
    )

    bed_trees = None
    if bed_file:
        bed_trees = load_bed_regions(bed_file)

    passed_breakpoints = []
    min_confidence = 0.7

    for i, bp in enumerate(breakpoints):
        bed_support_info = None
        if bed_trees:
            bed_support_info = get_support_read(bp['start_chr'], bp['start_loc'], bed_trees)

        tumor_coverage = min(ensure_int(bp['region_total_depths']['tumor'][1][0]),
                             ensure_int(bp['region_total_depths']['tumor'][1][1]))
        normal_coverage = min(ensure_int(bp['region_total_depths']['normal'][1][0]),
                              ensure_int(bp['region_total_depths']['normal'][1][1]))

        features = calculate_interpretable_features(bp, tumor_coverage, normal_coverage, platform)

        passed, reason = somatic_filter.filter_breakpoint(bp, bed_support_info)
        if not passed:
            continue

        confidence, confidence_factors = calculate_relative_confidence(features, reference_stats)

        if confidence >= min_confidence:
            bp['count'] = i
            bp['CONFIDENCE_SCORE'] = confidence
            bp['CONFIDENCE_FACTORS'] = confidence_factors
            bp['CLASS'] = 'PASSED_SOMATIC_FILTER'
            bp['FILTER_REASON'] = reason
            passed_breakpoints.append(bp)

    output_filtering_stats(passed_breakpoints, reference_stats, outdir)
    somatic_vcf_filename = outdir + ".only_long.somatic.vcf"
    write_somatic_vcf(passed_breakpoints, somatic_vcf_filename, ref_fasta_path)
    return len(passed_breakpoints), passed_breakpoints