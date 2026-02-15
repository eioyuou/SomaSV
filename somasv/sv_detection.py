import numpy as np
import pysam
from scipy import stats
from collections import defaultdict

from somasv.utils import validate_coordinates


def validate_coordinates(start, end):

    if start > end:
        return end, start
    return start, end

def detect_deletion_comprehensive(bam_file, chrom, start, end, global_cov,
                                  extend_bp=2000, bin_size=10,
                                  del_threshold=0.5, flanking_ratio=0.7):
    """
    Comprehensive deletion detection - original function preserved
    """

    print(f"Starting analysis for region: {chrom}:{start}-{end}")
    print(f"Global coverage: {global_cov:.2f}")
    print("-" * 50)

    # Step 1: Basic coverage analysis
    basic_results = calculate_basic_coverage(bam_file, chrom, start, end, global_cov,
                                             extend_bp, bin_size)
    if not basic_results:
        return None

    # Step 2: Statistical analysis
    stats_results = perform_statistical_analysis(basic_results)

    # Step 3: Gradient analysis
    gradient_results = analyze_coverage_gradient(basic_results)

    # Step 4: Anomaly detection
    anomaly_results = detect_coverage_anomalies(basic_results, global_cov)

    # Step 5: Comprehensive assessment
    final_assessment = comprehensive_assessment(
        basic_results, stats_results, gradient_results, anomaly_results,
        del_threshold, flanking_ratio
    )

    # Combine all results
    comprehensive_results = {
        'region_info': f"{chrom}:{start}-{end}",
        'basic_stats': basic_results,
        'statistical_analysis': stats_results,
        'gradient_analysis': gradient_results,
        'anomaly_detection': anomaly_results,
        'final_assessment': final_assessment
    }

    return comprehensive_results


def calculate_basic_coverage(bam_file, chrom, start, end, global_cov, extend_bp, bin_size):
    """Basic coverage calculation - original function preserved"""

    # Extend region around the target
    extended_start = max(1, start - extend_bp)
    extended_end = end + extend_bp

    # Open BAM file
    try:
        bamfile = pysam.AlignmentFile(bam_file, "rb")
    except Exception as e:
        print(f"Failed to open BAM file: {e}")
        return None

    try:
        coverage_arrays = bamfile.count_coverage(chrom, extended_start - 1, extended_end)
        total_coverage = np.sum(coverage_arrays, axis=0)
    except Exception as e:
        print(f"Failed to retrieve coverage: {e}")
        bamfile.close()
        return None

    bamfile.close()

    bin_coverages = []
    bin_positions = []

    for i in range(0, len(total_coverage), bin_size):
        bin_end = min(i + bin_size, len(total_coverage))
        bin_cov = np.mean(total_coverage[i:bin_end])
        bin_pos = extended_start + i + (bin_end - i) // 2

        bin_coverages.append(bin_cov)
        bin_positions.append(bin_pos)

    target_bins = []
    upstream_bins = []
    downstream_bins = []

    for i, pos in enumerate(bin_positions):
        if start <= pos <= end:
            target_bins.append(i)
        elif pos < start:
            upstream_bins.append(i)
        elif pos > end:
            downstream_bins.append(i)

    target_coverages = [bin_coverages[i] for i in target_bins] if target_bins else [0]
    upstream_coverages = [bin_coverages[i] for i in upstream_bins] if upstream_bins else []
    downstream_coverages = [bin_coverages[i] for i in downstream_bins] if downstream_bins else []
    flanking_coverages = upstream_coverages + downstream_coverages

    return {
        'positions': bin_positions,
        'coverages': bin_coverages,
        'target_bins': target_bins,
        'upstream_bins': upstream_bins,
        'downstream_bins': downstream_bins,
        'target_coverages': target_coverages,
        'flanking_coverages': flanking_coverages,
        'target_mean': np.mean(target_coverages),
        'target_median': np.median(target_coverages),
        'target_std': np.std(target_coverages),
        'flanking_mean': np.mean(flanking_coverages) if flanking_coverages else global_cov,
        'flanking_std': np.std(flanking_coverages) if len(flanking_coverages) > 1 else 0,
        'global_cov': global_cov
    }


def perform_statistical_analysis(basic_results):
    """统计分析 - 你的原函数保持不变"""
    target_coverages = basic_results['target_coverages']
    flanking_coverages = basic_results['flanking_coverages']
    global_cov = basic_results['global_cov']
    target_mean = basic_results['target_mean']
    flanking_mean = basic_results['flanking_mean']


    global_ratio = target_mean / global_cov if global_cov > 0 else 0
    flanking_ratio = target_mean / flanking_mean if flanking_mean > 0 else 0



    if len(target_coverages) > 1:
        t_stat_global, p_val_global = stats.ttest_1samp(target_coverages, global_cov)
    else:
        t_stat_global, p_val_global = 0, 1


    if len(target_coverages) > 1 and len(flanking_coverages) > 1:
        t_stat_flanking, p_val_flanking = stats.ttest_ind(target_coverages, flanking_coverages)
    else:
        t_stat_flanking, p_val_flanking = 0, 1


    all_coverages = np.array(basic_results['coverages'])
    z_scores = stats.zscore(all_coverages) if len(all_coverages) > 1 else np.zeros(len(all_coverages))
    target_z_scores = [z_scores[i] for i in basic_results['target_bins']] if basic_results['target_bins'] else [0]
    mean_z_score = np.mean(target_z_scores)

    # 4. 变异系数
    cv_target = basic_results['target_std'] / target_mean if target_mean > 0 else 0
    cv_flanking = basic_results['flanking_std'] / flanking_mean if flanking_mean > 0 else 0

    return {
        'global_ratio': global_ratio,
        'flanking_ratio': flanking_ratio,
        'p_value_vs_global': p_val_global,
        'p_value_vs_flanking': p_val_flanking,
        't_stat_global': t_stat_global,
        't_stat_flanking': t_stat_flanking,
        'mean_z_score': mean_z_score,
        'target_z_scores': target_z_scores,
        'cv_target': cv_target,
        'cv_flanking': cv_flanking
    }


def analyze_coverage_gradient(basic_results):

    positions = basic_results['positions']
    coverages = basic_results['coverages']
    target_bins = basic_results['target_bins']

    if not target_bins or len(target_bins) < 3:
        return {'has_gradient': False, 'gradient_type': 'insufficient_data'}


    window_size = 5
    pre_start = max(0, target_bins[0] - window_size)
    pre_end = target_bins[0]
    post_start = target_bins[-1] + 1
    post_end = min(len(coverages), target_bins[-1] + 1 + window_size)

    pre_cov = np.mean(coverages[pre_start:pre_end]) if pre_start < pre_end else basic_results['target_mean']
    target_cov = basic_results['target_mean']
    post_cov = np.mean(coverages[post_start:post_end]) if post_start < post_end else basic_results['target_mean']


    drop_ratio = target_cov / pre_cov if pre_cov > 0 else 1
    recovery_ratio = post_cov / target_cov if target_cov > 0 else 1

    # 判断梯度类型
    if drop_ratio < 0.7 and recovery_ratio > 1.3:
        gradient_type = 'sharp_drop_recovery'
        has_gradient = True
    elif drop_ratio < 0.8 and recovery_ratio > 1.2:
        gradient_type = 'moderate_drop_recovery'
        has_gradient = True
    elif 0.8 <= drop_ratio <= 1.2 and 0.8 <= recovery_ratio <= 1.2:
        gradient_type = 'smooth_transition'
        has_gradient = False
    else:
        gradient_type = 'irregular'
        has_gradient = False

    return {
        'has_gradient': has_gradient,
        'gradient_type': gradient_type,
        'pre_coverage': pre_cov,
        'target_coverage': target_cov,
        'post_coverage': post_cov,
        'drop_ratio': drop_ratio,
        'recovery_ratio': recovery_ratio
    }


def detect_coverage_anomalies(basic_results, global_cov):
    """Detect coverage anomalies - original deletion-oriented version"""

    target_mean = basic_results['target_mean']
    flanking_mean = basic_results['flanking_mean']

    anomalies = []

    # 1. Check if the flanking region is abnormally high
    flanking_vs_global = flanking_mean / global_cov
    if flanking_vs_global > 1.5:
        anomalies.append({
            'type': 'high_flanking_coverage',
            'description': f'Flanking region coverage is abnormally high ({flanking_vs_global:.2f}x of global)',
            'severity': 'high' if flanking_vs_global > 2.0 else 'medium'
        })

    # 2. Check target region coverage relative to global
    target_vs_global = target_mean / global_cov
    if target_vs_global > 1.2:
        anomalies.append({
            'type': 'high_target_coverage',
            'description': f'Target region coverage is higher than global ({target_vs_global:.2f}x)',
            'severity': 'low'
        })
    elif 0.8 <= target_vs_global <= 1.2:
        anomalies.append({
            'type': 'normal_target_coverage',
            'description': f'Target region coverage is close to global average ({target_vs_global:.2f}x)',
            'severity': 'info'
        })

    # 3. Identify lowest coverage continuous sub-region
    min_coverage_region = find_minimum_coverage_region(basic_results)

    return {
        'anomalies': anomalies,
        'flanking_vs_global_ratio': flanking_vs_global,
        'target_vs_global_ratio': target_vs_global,
        'min_coverage_region': min_coverage_region
    }


def find_minimum_coverage_region(basic_results, window_size=5):

    target_bins = basic_results['target_bins']
    coverages = basic_results['coverages']
    positions = basic_results['positions']

    if not target_bins or len(target_bins) < window_size:
        return None

    min_mean = float('inf')
    min_region = None

    for i in range(len(target_bins) - window_size + 1):
        window_bins = target_bins[i:i + window_size]
        window_coverages = [coverages[j] for j in window_bins]
        window_mean = np.mean(window_coverages)

        if window_mean < min_mean:
            min_mean = window_mean
            start_pos = positions[window_bins[0]]
            end_pos = positions[window_bins[-1]]
            min_region = {
                'start': start_pos,
                'end': end_pos,
                'mean_coverage': window_mean,
                'bins': window_bins
            }

    return min_region


def comprehensive_assessment(basic_results, stats_results, gradient_results,
                             anomaly_results, del_threshold, flanking_ratio):
    """Comprehensive assessment - original deletion detection logic preserved"""

    # Collect evidence
    strong_evidence = []
    medium_evidence = []
    weak_evidence = []
    counter_evidence = []

    target_mean = basic_results['target_mean']
    global_cov = basic_results['global_cov']

    # Strong evidence (+3 each)
    if stats_results['global_ratio'] < del_threshold:
        strong_evidence.append(
            f"Target coverage is lower than global by more than {del_threshold * 100:.0f}% (observed: {stats_results['global_ratio']:.1%})")

    if (stats_results['p_value_vs_global'] < 0.01 and
            target_mean < global_cov and
            stats_results['global_ratio'] < 0.8):
        strong_evidence.append(
            f"Highly significant difference from global coverage (p={stats_results['p_value_vs_global']:.2e})")

    if gradient_results['gradient_type'] == 'sharp_drop_recovery':
        strong_evidence.append("Sharp coverage drop and recovery pattern observed")

    # Medium evidence (+2 each)
    if (stats_results['flanking_ratio'] < flanking_ratio and
            stats_results['p_value_vs_flanking'] < 0.05 and
            not any(a['type'] == 'high_flanking_coverage' for a in anomaly_results['anomalies'])):
        medium_evidence.append(
            f"Significant drop relative to flanking regions (p={stats_results['p_value_vs_flanking']:.2e})")

    if stats_results['mean_z_score'] < -2:
        medium_evidence.append(f"Mean Z-score is significantly low ({stats_results['mean_z_score']:.2f})")

    if gradient_results['gradient_type'] == 'moderate_drop_recovery':
        medium_evidence.append("Moderate coverage drop and recovery pattern observed")

    # Weak evidence (+1 each)
    if del_threshold < stats_results['global_ratio'] < 0.8:
        weak_evidence.append(
            f"Target region coverage is slightly decreased relative to global ({stats_results['global_ratio']:.1%})")

    # Counter-evidence (-2 each)
    for anomaly in anomaly_results['anomalies']:
        if anomaly['type'] == 'high_flanking_coverage' and anomaly['severity'] == 'high':
            counter_evidence.append("Flanking coverage is abnormally high, may affect interpretation")
        elif anomaly['type'] == 'normal_target_coverage':
            counter_evidence.append("Target region coverage is close to global average")

    if gradient_results['gradient_type'] == 'smooth_transition':
        counter_evidence.append("Coverage changes smoothly, not consistent with deletion pattern")

    # Calculate confidence score
    confidence_score = (len(strong_evidence) * 3 +
                        len(medium_evidence) * 2 +
                        len(weak_evidence) * 1 -
                        len(counter_evidence) * 2)

    # Final interpretation
    if confidence_score >= 6 and len(strong_evidence) >= 1:
        conclusion = "Highly likely deletion variant"
        confidence_level = "high"
        is_deletion = True
    elif confidence_score >= 4 and len(strong_evidence) >= 1:
        conclusion = "Possible deletion variant"
        confidence_level = "medium"
        is_deletion = True
    elif confidence_score >= 2:
        conclusion = "Suspected deletion variant, validation recommended"
        confidence_level = "low"
        is_deletion = True
    else:
        conclusion = "More likely technical artifact or normal regional variation"
        confidence_level = "not_supported"
        is_deletion = False

    return {
        'is_deletion': is_deletion,
        'conclusion': conclusion,
        'confidence_level': confidence_level,
        'confidence_score': confidence_score,
        'max_possible_score': 12,  # Assume up to 4 strong evidence items
        'strong_evidence': strong_evidence,
        'medium_evidence': medium_evidence,
        'weak_evidence': weak_evidence,
        'counter_evidence': counter_evidence
    }




def get_reference_kmers(ref_fasta, chrom, start, end, k=31):

    try:
        seq = ref_fasta.fetch(chrom, start, end).upper()
        if len(seq) < k:
            return set()
        kmers = set(seq[i:i + k] for i in range(len(seq) - k + 1) if 'N' not in seq[i:i + k])
        return kmers
    except Exception as e:
        print(f"Error getting reference kmers: {e}")
        return set()


def count_kmers_from_reads(candidate_reads, k=31):

    kmer_counts = defaultdict(int)
    for read in candidate_reads:
        if read is None:
            continue
        seq = read.upper()
        if len(seq) < k:
            continue
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            if 'N' not in kmer:
                kmer_counts[kmer] += 1
    return kmer_counts


def match_kmers(reference_kmers, read_kmer_counts):

    matched = []
    for kmer in reference_kmers:
        count = read_kmer_counts.get(kmer, 0)
        matched.append((kmer, count))
    return matched


def is_kmer_match(ref, chrom, start, end, candidates):

    if not candidates:
        return False

    ref_kmers = get_reference_kmers(ref, chrom, start, end)
    read_kmers = count_kmers_from_reads(candidates, k=31)
    results = match_kmers(ref_kmers, read_kmers)


    well_supported_kmers = sum(1 for _, count in results if count >= 3)  # 至少3个reads支持
    total_kmers = len(results)

    if total_kmers == 0:
        return False

    fraction_supported = well_supported_kmers / total_kmers


    if fraction_supported > 0.85:
        return True  # germline
    else:
        return False  # somatic


def collect_kmer_candidates(bam_normal, chrom, start, end):

    candidates_kmer = []


    extended_start = max(0, start - 200)
    extended_end = end + 200

    for read in bam_normal.fetch(chrom, extended_start, extended_end):
        if read.cigartuples is None:
            continue

        # 收集多种类型的reads
        should_include = False


        if any(op == 4 and length >= 10 for op, length in read.cigartuples):
            should_include = True


        if any(op in [1, 2] and length >= 5 for op, length in read.cigartuples):
            should_include = True


        if read.mapping_quality < 20:
            should_include = True

        if should_include and read.query_sequence:
            candidates_kmer.append(read.query_sequence)

    return candidates_kmer


def evaluate_germline_evidence(bam_normal, ref_file, sv_candidate, count, margin=500):

    chrom_start = sv_candidate['start_chr']
    start_loc = sv_candidate['start_loc']
    end_loc = sv_candidate['end_loc']


    candidates_kmer = collect_kmer_candidates(bam_normal, chrom_start, start_loc, end_loc)


    evidence_score = 0
    evidence_details = []

    if len(candidates_kmer) >= 2:
        is_kmer_match_result = is_kmer_match(ref_file, chrom_start, start_loc, end_loc, candidates_kmer)
        if is_kmer_match_result:
            evidence_score += 2
            evidence_details.append("K-mer: Low coverage detected")

    return evidence_score >= 2, evidence_score, evidence_details, len(candidates_kmer)



def detect_duplication_comprehensive(bam_file, chrom, start, end, global_cov,
                                     extend_bp=2000, bin_size=10,
                                     dup_threshold=1.5, flanking_ratio=1.3):
    """
    Comprehensive duplication detection - based on deletion logic with reversed thresholds.
    """

    print(f"Starting duplication analysis: {chrom}:{start}-{end}")
    print(f"Global coverage: {global_cov:.2f}")
    print("-" * 50)

    # Step 1: Basic coverage analysis (reusing your function)
    basic_results = calculate_basic_coverage(bam_file, chrom, start, end, global_cov,
                                             extend_bp, bin_size)
    if not basic_results:
        return None

    # Step 2: Statistical analysis (reusing your function)
    stats_results = perform_statistical_analysis(basic_results)

    # Step 3: Gradient analysis (adjusted for duplication patterns)
    gradient_results = analyze_coverage_gradient_dup(basic_results)

    # Step 4: Anomaly detection (adjusted for duplication logic)
    anomaly_results = detect_coverage_anomalies_dup(basic_results, global_cov)

    # Step 5: Comprehensive assessment (specific to duplication)
    final_assessment = comprehensive_assessment_dup(
        basic_results, stats_results, gradient_results, anomaly_results,
        dup_threshold, flanking_ratio
    )

    # Combine all results
    comprehensive_results = {
        'region_info': f"{chrom}:{start}-{end}",
        'basic_stats': basic_results,
        'statistical_analysis': stats_results,
        'gradient_analysis': gradient_results,
        'anomaly_detection': anomaly_results,
        'final_assessment': final_assessment
    }

    return comprehensive_results


def analyze_coverage_gradient_dup(basic_results):
    """Duplication coverage gradient analysis - opposite of deletion"""
    positions = basic_results['positions']
    coverages = basic_results['coverages']
    target_bins = basic_results['target_bins']

    if not target_bins or len(target_bins) < 3:
        return {'has_gradient': False, 'gradient_type': 'insufficient_data'}

    # Analyze coverage before and after the target region
    window_size = 5
    pre_start = max(0, target_bins[0] - window_size)
    pre_end = target_bins[0]
    post_start = target_bins[-1] + 1
    post_end = min(len(coverages), target_bins[-1] + 1 + window_size)

    pre_cov = np.mean(coverages[pre_start:pre_end]) if pre_start < pre_end else basic_results['target_mean']
    target_cov = basic_results['target_mean']
    post_cov = np.mean(coverages[post_start:post_end]) if post_start < post_end else basic_results['target_mean']

    # Calculate elevation ratio (reverse of deletion logic)
    elevation_ratio_pre = target_cov / pre_cov if pre_cov > 0 else 1
    elevation_ratio_post = target_cov / post_cov if post_cov > 0 else 1

    # Determine gradient type (duplication should show elevated pattern)
    if elevation_ratio_pre > 1.3 and elevation_ratio_post > 1.3:
        gradient_type = 'sharp_elevation_pattern'  # Sharp increase, typical duplication
        has_gradient = True
    elif elevation_ratio_pre > 1.2 and elevation_ratio_post > 1.2:
        gradient_type = 'moderate_elevation_pattern'  # Moderate increase
        has_gradient = True
    elif 0.8 <= elevation_ratio_pre <= 1.2 and 0.8 <= elevation_ratio_post <= 1.2:
        gradient_type = 'smooth_transition'  # Smooth transition
        has_gradient = False
    else:
        gradient_type = 'irregular'  # Irregular change
        has_gradient = False

    return {
        'has_gradient': has_gradient,
        'gradient_type': gradient_type,
        'pre_coverage': pre_cov,
        'target_coverage': target_cov,
        'post_coverage': post_cov,
        'elevation_ratio_pre': elevation_ratio_pre,
        'elevation_ratio_post': elevation_ratio_post
    }


def detect_coverage_anomalies_dup(basic_results, global_cov):
    """Duplication anomaly detection – opposite logic to deletions"""
    target_mean = basic_results['target_mean']
    flanking_mean = basic_results['flanking_mean']

    anomalies = []

    # 1. Check target region coverage relative to global (duplication should show higher coverage)
    target_vs_global = target_mean / global_cov
    if target_vs_global > 2.0:
        anomalies.append({
            'type': 'very_high_target_coverage',
            'description': f'Target region coverage is significantly higher than global ({target_vs_global:.2f}x)',
            'severity': 'high'
        })
    elif target_vs_global > 1.5:
        anomalies.append({
            'type': 'high_target_coverage',
            'description': f'Target region coverage is elevated compared to global ({target_vs_global:.2f}x)',
            'severity': 'medium'
        })
    elif 0.8 <= target_vs_global <= 1.2:
        anomalies.append({
            'type': 'normal_target_coverage',
            'description': f'Target region coverage is close to global expectation ({target_vs_global:.2f}x)',
            'severity': 'info'
        })

    # 2. Check if flanking regions have abnormally low coverage (could affect interpretation)
    flanking_vs_global = flanking_mean / global_cov
    if flanking_vs_global < 0.5:
        anomalies.append({
            'type': 'low_flanking_coverage',
            'description': f'Flanking region coverage is abnormally low ({flanking_vs_global:.2f}x of global)',
            'severity': 'medium'
        })

    return {
        'anomalies': anomalies,
        'flanking_vs_global_ratio': flanking_vs_global,
        'target_vs_global_ratio': target_vs_global
    }


def comprehensive_assessment_dup(basic_results, stats_results, gradient_results,
                                 anomaly_results, dup_threshold, flanking_ratio):
    """Comprehensive duplication assessment - opposite to deletions"""
    # Collect evidence
    strong_evidence = []
    medium_evidence = []
    weak_evidence = []
    counter_evidence = []

    target_mean = basic_results['target_mean']
    global_cov = basic_results['global_cov']

    # Strong evidence (+3 each) - duplication should show high coverage
    if stats_results['global_ratio'] > dup_threshold:
        strong_evidence.append(f"Target region coverage is higher than global by a factor of {dup_threshold} (observed: {stats_results['global_ratio']:.2f}x)")

    if (stats_results['p_value_vs_global'] < 0.01 and
            target_mean > global_cov and
            stats_results['global_ratio'] > 1.2):
        strong_evidence.append(f"Highly significant difference from global coverage (p={stats_results['p_value_vs_global']:.2e})")

    if gradient_results['gradient_type'] == 'sharp_elevation_pattern':
        strong_evidence.append("Sharp elevation pattern in coverage observed")

    # Medium evidence (+2 each)
    if (stats_results['flanking_ratio'] > flanking_ratio and
            stats_results['p_value_vs_flanking'] < 0.05 and
            not any(a['type'] == 'low_flanking_coverage' for a in anomaly_results['anomalies'])):
        medium_evidence.append(f"Significant increase relative to flanking regions (p={stats_results['p_value_vs_flanking']:.2e})")

    if stats_results['mean_z_score'] > 2:
        medium_evidence.append(f"Mean Z-score is noticeably high ({stats_results['mean_z_score']:.2f})")

    if gradient_results['gradient_type'] == 'moderate_elevation_pattern':
        medium_evidence.append("Moderate elevation pattern in coverage observed")

    # Weak evidence (+1 each)
    if 1.2 < stats_results['global_ratio'] < dup_threshold:
        weak_evidence.append(f"Target region coverage is slightly higher than global ({stats_results['global_ratio']:.2f}x)")

    # Counter-evidence (-2 each)
    for anomaly in anomaly_results['anomalies']:
        if anomaly['type'] == 'low_flanking_coverage' and anomaly['severity'] == 'medium':
            counter_evidence.append("Flanking region has abnormally low coverage, may confound interpretation")
        elif anomaly['type'] == 'normal_target_coverage':
            counter_evidence.append("Target region coverage is similar to normal global level")

    if gradient_results['gradient_type'] == 'smooth_transition':
        counter_evidence.append("Coverage transitions smoothly; not fitting duplication pattern")

    # Compute confidence score
    confidence_score = (len(strong_evidence) * 3 +
                        len(medium_evidence) * 2 +
                        len(weak_evidence) * 1 -
                        len(counter_evidence) * 2)

    # Final decision
    if confidence_score >= 6 and len(strong_evidence) >= 1:
        conclusion = "Highly likely duplication event"
        confidence_level = "high"
        is_duplication = True
    elif confidence_score >= 4 and len(strong_evidence) >= 1:
        conclusion = "Possible duplication event"
        confidence_level = "medium"
        is_duplication = True
    elif confidence_score >= 2:
        conclusion = "Suspected duplication, further validation recommended"
        confidence_level = "low"
        is_duplication = True
    else:
        conclusion = "More likely technical artifact or regional variability"
        confidence_level = "not_supported"
        is_duplication = False

    return {
        'is_duplication': is_duplication,
        'conclusion': conclusion,
        'confidence_level': confidence_level,
        'confidence_score': confidence_score,
        'max_possible_score': 12,
        'strong_evidence': strong_evidence,
        'medium_evidence': medium_evidence,
        'weak_evidence': weak_evidence,
        'counter_evidence': counter_evidence
    }


def is_kmer_match_dup(ref, chrom, start, end, candidates):

    if not candidates:
        return False

    ref_kmers = get_reference_kmers(ref, chrom, start, end)
    read_kmers = count_kmers_from_reads(candidates, k=31)


    high_frequency_kmers = 0
    total_kmers = len(ref_kmers)

    if total_kmers == 0:
        return False

    # 估算预期频率
    expected_freq = len(candidates) / max(100, total_kmers)  # 粗略估计

    for kmer in ref_kmers:
        freq = read_kmers.get(kmer, 0)
        if freq > expected_freq * 1.5:
            high_frequency_kmers += 1

    fraction_high_freq = high_frequency_kmers / total_kmers

    print(f"DUP K-mer analysis: {high_frequency_kmers}/{total_kmers} high-frequency k-mers ({fraction_high_freq:.3f})")


    if fraction_high_freq > 0.3:
        return True  # germline
    else:
        return False  # somatic


def evaluate_germline_evidence_dup(bam_normal, ref_file, sv_candidate, count, margin=500):

    chrom_start = sv_candidate['start_chr']
    start_loc = sv_candidate['start_loc']
    end_loc = sv_candidate['end_loc']

    # 收集k-mer分析用的reads
    candidates_kmer = collect_kmer_candidates(bam_normal, chrom_start, start_loc, end_loc)

    # 计算证据得分
    evidence_score = 0
    evidence_details = []

    if len(candidates_kmer) >= 2:
        is_kmer_match_result = is_kmer_match_dup(ref_file, chrom_start, start_loc, end_loc, candidates_kmer)
        if is_kmer_match_result:
            evidence_score += 2
            evidence_details.append("K-mer: High frequency detected")

    return evidence_score >= 2, evidence_score, evidence_details, len(candidates_kmer)


def detect_inversion_validation(bam_file, chrom, start, end, min_mapq=20):
    """
    Improved inversion validation - added more evidence layers
    """
    print(f"Inversion validation: {chrom}:{start}-{end}")

    start, end = validate_coordinates(start, end)
    region_size = end - start

    # Core statistics
    total_pairs = 0
    same_orientation_pairs = 0
    opposite_orientation_pairs = 0

    # New: Split reads near breakpoints
    left_breakpoint_splits = []
    right_breakpoint_splits = []

    # New: Insert size analysis
    insert_sizes = []
    abnormal_insert_sizes = 0

    # New: Read strand direction inside inversion
    forward_reads = 0
    reverse_reads = 0

    # Improved search window
    search_start = start + min(100, region_size // 10)
    search_end = end - min(100, region_size // 10)

    if search_end <= search_start:
        search_start = max(0, start - 500)
        search_end = end + 500

    try:
        bamfile = pysam.AlignmentFile(bam_file, "rb")

        # 1. Analyze paired-end reads inside inversion region
        for read in bamfile.fetch(chrom, search_start, search_end):
            if (read.is_paired and not read.is_unmapped and not read.mate_is_unmapped and
                    read.reference_name == read.next_reference_name and
                    read.mapping_quality >= min_mapq and
                    not read.is_secondary and not read.is_supplementary):

                mate_pos = read.next_reference_start
                if not (search_start <= mate_pos <= search_end):
                    continue

                total_pairs += 1

                # Same or opposite strand
                if read.is_reverse == read.mate_is_reverse:
                    same_orientation_pairs += 1
                else:
                    opposite_orientation_pairs += 1

                # Insert size check
                if hasattr(read, 'template_length') and read.template_length != 0:
                    insert_size = abs(read.template_length)
                    insert_sizes.append(insert_size)
                    if insert_size > 2000 or insert_size < 100:
                        abnormal_insert_sizes += 1

            # Strand bias stats
            if (not read.is_unmapped and read.mapping_quality >= min_mapq and
                    not read.is_secondary and not read.is_supplementary):
                if read.is_reverse:
                    reverse_reads += 1
                else:
                    forward_reads += 1

        # 2. Enhanced split read detection around breakpoints (soft-clips, indels)
        for breakpoint, split_list in [(start, left_breakpoint_splits), (end, right_breakpoint_splits)]:
            for read in bamfile.fetch(chrom, max(0, breakpoint - 200), breakpoint + 200):
                if read.mapping_quality >= min_mapq and read.cigartuples:

                    # Soft clipping and complex CIGAR operations
                    soft_clips = []
                    has_complex_cigar = False

                    for op, length in read.cigartuples:
                        if op == 4 and length >= 15:
                            soft_clips.append(length)
                        elif op in [1, 2] and length >= 10:
                            has_complex_cigar = True

                    if soft_clips and abs(read.reference_start - breakpoint) <= 100:
                        split_info = {
                            'pos': read.reference_start,
                            'clip_lengths': soft_clips,
                            'distance_to_breakpoint': abs(read.reference_start - breakpoint),
                            'mapq': read.mapping_quality,
                            'has_complex_cigar': has_complex_cigar
                        }
                        split_list.append(split_info)

        # 3. Detect supplementary alignments that may indicate inversion
        supplementary_evidence = 0
        for read in bamfile.fetch(chrom, start - 1000, end + 1000):
            if read.is_supplementary and read.mapping_quality >= min_mapq:
                if read.has_tag('SA'):
                    sa_tag = read.get_tag('SA')
                    if chrom in sa_tag and ('+-' in sa_tag or '-+' in sa_tag):
                        supplementary_evidence += 1

        bamfile.close()

    except Exception as e:
        print(f"Inversion validation failed: {e}")
        return {
            'is_real_inversion': False,
            'conclusion': f'Validation failed: {e}',
            'confidence_level': 'error'
        }

    if total_pairs == 0:
        return {
            'is_real_inversion': False,
            'conclusion': 'No paired-end reads found in region',
            'confidence_level': 'insufficient_data'
        }

    # Calculate ratios and metrics
    same_orientation_ratio = same_orientation_pairs / total_pairs

    total_single_reads = forward_reads + reverse_reads
    strand_bias = abs(forward_reads - reverse_reads) / total_single_reads if total_single_reads > 0 else 0

    abnormal_insert_ratio = abnormal_insert_sizes / len(insert_sizes) if insert_sizes else 0

    print(f"Inversion analysis summary:")
    print(f"  Total pairs: {total_pairs}")
    print(f"  Same-orientation pairs: {same_orientation_pairs} ({same_orientation_ratio:.2%})")
    print(f"  Opposite-orientation pairs: {opposite_orientation_pairs}")
    print(f"  Split reads at left breakpoint: {len(left_breakpoint_splits)}")
    print(f"  Split reads at right breakpoint: {len(right_breakpoint_splits)}")
    print(f"  Abnormal insert sizes: {abnormal_insert_sizes}/{len(insert_sizes)} ({abnormal_insert_ratio:.2%})")
    print(f"  Strand bias: {strand_bias:.2%}")
    print(f"  Supplementary alignment evidence: {supplementary_evidence}")

    # Scoring evidences
    evidence_score = 0
    evidence_details = []

    if same_orientation_ratio >= 0.7:
        evidence_score += 5
        evidence_details.append(f"High ratio of same-orientation pairs ({same_orientation_ratio:.1%})")
    elif same_orientation_ratio >= 0.5:
        evidence_score += 4
        evidence_details.append(f"Moderate same-orientation ratio ({same_orientation_ratio:.1%})")
    elif same_orientation_ratio >= 0.3:
        evidence_score += 2
        evidence_details.append(f"Some same-orientation support ({same_orientation_ratio:.1%})")
    elif same_orientation_ratio <= 0.1:
        evidence_score -= 3
        evidence_details.append(f"Low same-orientation support ({same_orientation_ratio:.1%})")

    if same_orientation_pairs >= 8:
        evidence_score += 3
        evidence_details.append(f"Sufficient number of same-orientation pairs ({same_orientation_pairs})")
    elif same_orientation_pairs >= 3:
        evidence_score += 2
        evidence_details.append(f"Some same-orientation pairs ({same_orientation_pairs})")
    elif same_orientation_pairs >= 1:
        evidence_score += 1
        evidence_details.append(f"A few same-orientation pairs ({same_orientation_pairs})")

    total_splits = len(left_breakpoint_splits) + len(right_breakpoint_splits)
    if total_splits >= 4:
        evidence_score += 3
        evidence_details.append(f"Multiple breakpoint split reads ({total_splits})")
    elif total_splits >= 2:
        evidence_score += 2
        evidence_details.append(f"A few breakpoint split reads ({total_splits})")
    elif total_splits >= 1:
        evidence_score += 1
        evidence_details.append(f"Rare split read support ({total_splits})")

    if abnormal_insert_ratio >= 0.5:
        evidence_score += 3
        evidence_details.append(f"High abnormal insert size ratio ({abnormal_insert_ratio:.1%})")
    elif abnormal_insert_ratio >= 0.3:
        evidence_score += 2
        evidence_details.append(f"Some abnormal insert sizes ({abnormal_insert_ratio:.1%})")

    if supplementary_evidence >= 3:
        evidence_score += 3
        evidence_details.append(f"Strong supplementary alignment evidence ({supplementary_evidence})")
    elif supplementary_evidence >= 1:
        evidence_score += 2
        evidence_details.append(f"Some supplementary alignment evidence ({supplementary_evidence})")

    if strand_bias >= 0.7:
        evidence_score += 2
        evidence_details.append(f"Significant strand bias ({strand_bias:.1%})")
    elif strand_bias >= 0.6:
        evidence_score += 1
        evidence_details.append(f"Moderate strand bias ({strand_bias:.1%})")

    # Final decision
    if evidence_score >= 8:
        is_real = True
        conclusion = f"Strong inversion signal: {'; '.join(evidence_details)}"
        confidence = "high"
    elif evidence_score >= 5:
        is_real = True
        conclusion = f"Moderate inversion signal: {'; '.join(evidence_details)}"
        confidence = "medium"
    elif evidence_score >= 2:
        is_real = True
        conclusion = f"Weak inversion signal: {'; '.join(evidence_details)}"
        confidence = "low"
    else:
        is_real = False
        conclusion = f"Insufficient inversion evidence: {'; '.join(evidence_details)}"
        confidence = "negative"

    return {
        'is_real_inversion': is_real,
        'conclusion': conclusion,
        'confidence_level': confidence,
        'evidence_score': evidence_score,
        'same_orientation_ratio': same_orientation_ratio,
        'same_orientation_pairs': same_orientation_pairs,
        'total_pairs': total_pairs,
        'split_reads_count': total_splits,
        'abnormal_insert_ratio': abnormal_insert_ratio,
        'supplementary_evidence': supplementary_evidence
    }

def detect_bnd_by_discordant_pairs(bam_file, chrom1, pos1, chrom2, pos2, window=5000, min_mapq=20):
    """
    Correct BND (translocation) detection logic
    """
    print(f"Detecting translocation: {chrom1}:{pos1} <-> {chrom2}:{pos2}")

    # Statistic counters
    discordant_pairs_1_to_2 = 0  # read on chrom1, mate on chrom2
    discordant_pairs_2_to_1 = 0  # read on chrom2, mate on chrom1
    split_reads_1_to_2 = 0       # split reads from chrom1 pointing to chrom2
    split_reads_2_to_1 = 0       # split reads from chrom2 pointing to chrom1

    supporting_reads_chrom1 = 0
    supporting_reads_chrom2 = 0

    # Orientation pattern tracking
    orientation_patterns = {'FF': 0, 'FR': 0, 'RF': 0, 'RR': 0}

    try:
        bamfile = pysam.AlignmentFile(bam_file, "rb")

        # === Analyze first breakpoint region (chrom1:pos1) ===
        print(f"Analyzing breakpoint 1: {chrom1}:{pos1 - window}-{pos1 + window}")

        for read in bamfile.fetch(chrom1, max(0, pos1 - window), pos1 + window):
            if read.is_unmapped or read.mapping_quality < min_mapq:
                continue
            if read.is_secondary or read.is_supplementary:
                continue

            supporting_reads_chrom1 += 1

            # 1. Check discordant pairs: read on chrom1, mate on chrom2
            if (read.is_paired and not read.mate_is_unmapped and
                    read.next_reference_name == chrom2):

                mate_distance_to_bp2 = abs(read.next_reference_start - pos2)
                if mate_distance_to_bp2 <= window:
                    discordant_pairs_1_to_2 += 1

                    # Record orientation pattern
                    read_orient = 'R' if read.is_reverse else 'F'
                    mate_orient = 'R' if read.mate_is_reverse else 'F'
                    pattern = read_orient + mate_orient
                    orientation_patterns[pattern] += 1

            # 2. Check for potential split reads
            if read.cigartuples and not read.is_paired:
                soft_clip_length = 0
                for op, length in read.cigartuples:
                    if op == 4 and length >= 20:
                        soft_clip_length = max(soft_clip_length, length)

                if soft_clip_length >= 20 and abs(read.reference_start - pos1) <= 100:
                    # Ideally check mapping of clipped sequence to chrom2:pos2
                    split_reads_1_to_2 += 1
                    print(f"  Found potential split read: {read.query_name} at {chrom1}:{read.reference_start}, soft-clipped {soft_clip_length}bp")

            # 3. Check supplementary alignments (SA tag)
            if read.has_tag('SA'):
                sa_tag = read.get_tag('SA')
                sa_parts = sa_tag.strip(';').split(';')
                for sa_part in sa_parts:
                    if not sa_part:
                        continue
                    try:
                        sa_fields = sa_part.split(',')
                        if len(sa_fields) >= 6:
                            sa_chrom = sa_fields[0]
                            sa_pos = int(sa_fields[1])
                            sa_mapq = int(sa_fields[4])

                            if (sa_chrom == chrom2 and
                                    abs(sa_pos - pos2) <= window and
                                    sa_mapq >= min_mapq):
                                split_reads_1_to_2 += 1
                                print(f"  Found split read: {read.query_name} {chrom1}:{read.reference_start} -> {chrom2}:{sa_pos}")
                    except (ValueError, IndexError):
                        continue

        # === Analyze second breakpoint region (chrom2:pos2) ===
        print(f"Analyzing breakpoint 2: {chrom2}:{pos2 - window}-{pos2 + window}")

        for read in bamfile.fetch(chrom2, max(0, pos2 - window), pos2 + window):
            if read.is_unmapped or read.mapping_quality < min_mapq:
                continue
            if read.is_secondary or read.is_supplementary:
                continue

            supporting_reads_chrom2 += 1

            # 1. Check discordant pairs: read on chrom2, mate on chrom1
            if (read.is_paired and not read.mate_is_unmapped and
                    read.next_reference_name == chrom1):

                mate_distance_to_bp1 = abs(read.next_reference_start - pos1)
                if mate_distance_to_bp1 <= window:
                    discordant_pairs_2_to_1 += 1
                    print(f"  Found discordant pair: {read.query_name} {chrom2}:{read.reference_start} -> {chrom1}:{read.next_reference_start}")

            # 2. Check split reads from chrom2 pointing to chrom1
            if read.has_tag('SA'):
                sa_tag = read.get_tag('SA')
                sa_parts = sa_tag.strip(';').split(';')
                for sa_part in sa_parts:
                    if not sa_part:
                        continue
                    try:
                        sa_fields = sa_part.split(',')
                        if len(sa_fields) >= 6:
                            sa_chrom = sa_fields[0]
                            sa_pos = int(sa_fields[1])
                            sa_mapq = int(sa_fields[4])

                            if (sa_chrom == chrom1 and
                                    abs(sa_pos - pos1) <= window and
                                    sa_mapq >= min_mapq):
                                split_reads_2_to_1 += 1
                                print(f"  Found split read: {read.query_name} {chrom2}:{read.reference_start} -> {chrom1}:{sa_pos}")
                    except (ValueError, IndexError):
                        continue

        bamfile.close()

    except Exception as e:
        print(f"Translocation detection failed: {e}")
        return {'is_bnd': False, 'conclusion': f'Analysis failed: {e}', 'confidence_level': 'error'}

    # Summary statistics
    total_discordant_pairs = discordant_pairs_1_to_2 + discordant_pairs_2_to_1
    total_split_reads = split_reads_1_to_2 + split_reads_2_to_1
    total_supporting_reads = supporting_reads_chrom1 + supporting_reads_chrom2

    # Analyze orientation patterns
    total_oriented_pairs = sum(orientation_patterns.values())
    dominant_pattern = max(orientation_patterns.items(), key=lambda x: x[1]) if total_oriented_pairs > 0 else ('None', 0)

    print(f"\nTranslocation analysis summary:")
    print(f"  Discordant pairs ({chrom1}->{chrom2}): {discordant_pairs_1_to_2}")
    print(f"  Discordant pairs ({chrom2}->{chrom1}): {discordant_pairs_2_to_1}")
    print(f"  Total discordant pairs: {total_discordant_pairs}")
    print(f"  Split reads ({chrom1}->{chrom2}): {split_reads_1_to_2}")
    print(f"  Split reads ({chrom2}->{chrom1}): {split_reads_2_to_1}")
    print(f"  Total split reads: {total_split_reads}")
    print(f"  Supporting reads: {total_supporting_reads}")
    print(f"  Dominant orientation pattern: {dominant_pattern[0]} ({dominant_pattern[1]} occurrences)")

    # Scoring evidence
    evidence_score = 0
    evidence_details = []

    # Discordant pairs
    if total_discordant_pairs >= 5:
        evidence_score += 6
        evidence_details.append(f"Strong discordant pair support ({total_discordant_pairs})")
    elif total_discordant_pairs >= 3:
        evidence_score += 4
        evidence_details.append(f"Moderate discordant pair support ({total_discordant_pairs})")
    elif total_discordant_pairs >= 1:
        evidence_score += 2
        evidence_details.append(f"Some discordant pair support ({total_discordant_pairs})")

    # Split reads
    if total_split_reads >= 3:
        evidence_score += 5
        evidence_details.append(f"Strong split read support ({total_split_reads})")
    elif total_split_reads >= 2:
        evidence_score += 3
        evidence_details.append(f"Moderate split read support ({total_split_reads})")
    elif total_split_reads >= 1:
        evidence_score += 2
        evidence_details.append(f"Some split read support ({total_split_reads})")

    # Bidirectional evidence
    if discordant_pairs_1_to_2 > 0 and discordant_pairs_2_to_1 > 0:
        evidence_score += 3
        evidence_details.append("Bidirectional discordant pair evidence")

    if split_reads_1_to_2 > 0 and split_reads_2_to_1 > 0:
        evidence_score += 3
        evidence_details.append("Bidirectional split read evidence")

    # Orientation consistency
    if dominant_pattern[1] >= 2 and total_oriented_pairs > 0:
        consistency = dominant_pattern[1] / total_oriented_pairs
        if consistency >= 0.7:
            evidence_score += 2
            evidence_details.append(f"Consistent orientation pattern ({dominant_pattern[0]}: {consistency:.1%})")

    # Regional support
    if total_supporting_reads >= 20:
        evidence_score += 1
        evidence_details.append(f"Adequate regional support ({total_supporting_reads} reads)")

    # Final evaluation
    if evidence_score >= 10:
        conclusion = f"Strong translocation signal: {'; '.join(evidence_details)}"
        confidence = 'high'
        is_bnd = True
    elif evidence_score >= 6:
        conclusion = f"Moderate translocation signal: {'; '.join(evidence_details)}"
        confidence = 'medium'
        is_bnd = True
    elif evidence_score >= 3:
        conclusion = f"Weak translocation signal: {'; '.join(evidence_details)}"
        confidence = 'low'
        is_bnd = True
    else:
        conclusion = f"No translocation evidence: total score {evidence_score}, {'; '.join(evidence_details) if evidence_details else 'no significant signal'}"
        confidence = 'not_supported'
        is_bnd = False

    return {
        'is_bnd': is_bnd,
        'conclusion': conclusion,
        'confidence_level': confidence,
        'evidence_score': evidence_score,
        'total_discordant_pairs': total_discordant_pairs,
        'total_split_reads': total_split_reads,
        'bidirectional_evidence': (discordant_pairs_1_to_2 > 0 and discordant_pairs_2_to_1 > 0) or (
                    split_reads_1_to_2 > 0 and split_reads_2_to_1 > 0)
    }


def detect_insertion_comprehensive(bam_file, chrom, pos, inserted_seq,
                                   window=500, min_mapq=20, min_clip_length=10):
    """
    Comprehensive insertion detection
    """
    print(f"Start analyzing insertion: {chrom}:{pos}")
    print(f"Inserted sequence: {inserted_seq[:50]}{'...' if len(inserted_seq) > 50 else ''} ({len(inserted_seq)}bp)")
    print("-" * 50)

    # Statistic variables
    total_reads = 0
    split_reads = 0
    supporting_soft_clips = 0

    # Soft-clipped sequences collection
    left_clips = []  # Soft clips on the left side of the breakpoint
    right_clips = []  # Soft clips on the right side of the breakpoint

    # k-mers from the inserted sequence (for validation)
    insertion_kmers = set()
    if len(inserted_seq) >= 21 and inserted_seq != 'N/A':
        k = min(21, len(inserted_seq) // 2)  # Dynamically adjust k-mer length
        insertion_kmers = {inserted_seq[i:i + k] for i in range(len(inserted_seq) - k + 1)
                           if 'N' not in inserted_seq[i:i + k]}

    try:
        if isinstance(bam_file, pysam.AlignmentFile):
            print("Error: Expected a file path, got an AlignmentFile object")
            return {'is_insertion': False, 'conclusion': 'Invalid parameter type', 'confidence_level': 'error'}

        try:
            # Open BAM file
            bamfile = pysam.AlignmentFile(bam_file, "rb")
        except Exception as e:
            return {'is_insertion': False, 'conclusion': f'Failed to open BAM file: {e}', 'confidence_level': 'error'}

        # Search for reads near the insertion site
        for read in bamfile.fetch(chrom, max(0, pos - window), pos + window):
            if (read.is_unmapped or read.mapping_quality < min_mapq or
                    read.is_secondary or read.is_supplementary):
                continue

            total_reads += 1

            # Check for soft clipping in CIGAR
            if read.cigartuples:
                read_start = read.reference_start
                read_end = read.reference_end

                # Check if the read is near the insertion site
                if abs(read_start - pos) <= 100 or abs(read_end - pos) <= 100:

                    # Analyze CIGAR string
                    cigar_ops = read.cigartuples
                    query_seq = read.query_sequence

                    if query_seq is None:
                        continue

                    # Check left soft clip (5' soft clip)
                    if cigar_ops[0][0] == 4 and cigar_ops[0][1] >= min_clip_length:
                        clip_length = cigar_ops[0][1]
                        clip_seq = query_seq[:clip_length]

                        # If read is to the right of the insertion site, left soft clip may contain the insertion
                        if read_start >= pos - 50:
                            left_clips.append({
                                'seq': clip_seq,
                                'length': clip_length,
                                'read_pos': read_start,
                                'distance_to_breakpoint': abs(read_start - pos)
                            })
                            split_reads += 1

                    # Check right soft clip (3' soft clip)
                    if cigar_ops[-1][0] == 4 and cigar_ops[-1][1] >= min_clip_length:
                        clip_length = cigar_ops[-1][1]
                        clip_seq = query_seq[-clip_length:]

                        # If read is to the left of the insertion site, right soft clip may contain the insertion
                        if read_end <= pos + 50:
                            right_clips.append({
                                'seq': clip_seq,
                                'length': clip_length,
                                'read_pos': read_end,
                                'distance_to_breakpoint': abs(read_end - pos)
                            })
                            split_reads += 1

        bamfile.close()

    except Exception as e:
        print(f"Insertion detection failed: {e}")
        return None

    # Analyze soft-clipped sequences
    insertion_evidence_score = 0
    evidence_details = []

    print(f"Detected soft-clipped sequences:")
    print(f"  Left soft clips: {len(left_clips)}")
    print(f"  Right soft clips: {len(right_clips)}")
    print(f"  Total split reads: {split_reads}")

    # K-mer matching analysis
    if insertion_kmers and (left_clips or right_clips):
        matched_clips = 0
        total_clips = len(left_clips) + len(right_clips)

        # Check left soft clips
        for clip in left_clips:
            if contains_insertion_kmers(clip['seq'], insertion_kmers):
                matched_clips += 1
                print(f"    Left soft clip matches insertion sequence: {clip['seq'][:30]}...")

        # Check right soft clips
        for clip in right_clips:
            if contains_insertion_kmers(clip['seq'], insertion_kmers):
                matched_clips += 1
                print(f"    Right soft clip matches insertion sequence: {clip['seq'][:30]}...")

        if total_clips > 0:
            match_ratio = matched_clips / total_clips
            print(f"  K-mer match ratio: {matched_clips}/{total_clips} ({match_ratio:.2%})")

            if match_ratio >= 0.5:
                insertion_evidence_score += 4
                evidence_details.append(f"High proportion of soft clips matched insertion sequence ({match_ratio:.1%})")
            elif match_ratio >= 0.3:
                insertion_evidence_score += 3
                evidence_details.append(f"Moderate proportion of soft clips matched insertion sequence ({match_ratio:.1%})")
            elif match_ratio >= 0.1:
                insertion_evidence_score += 2
                evidence_details.append(f"Low proportion of soft clips matched insertion sequence ({match_ratio:.1%})")

    # Split reads count evidence
    if split_reads >= 5:
        insertion_evidence_score += 3
        evidence_details.append(f"Multiple split reads support detected ({split_reads})")
    elif split_reads >= 2:
        insertion_evidence_score += 2
        evidence_details.append(f"A few split reads support detected ({split_reads})")
    elif split_reads >= 1:
        insertion_evidence_score += 1
        evidence_details.append(f"Limited split reads support ({split_reads})")

    # Bidirectional soft-clip support
    if len(left_clips) > 0 and len(right_clips) > 0:
        insertion_evidence_score += 2
        evidence_details.append("Bidirectional soft clip evidence")

    # Consistency between soft clip length and inserted sequence length
    if inserted_seq != 'N/A' and len(inserted_seq) > 10:
        expected_length = len(inserted_seq)
        length_consistent_clips = 0

        for clip in left_clips + right_clips:
            if abs(clip['length'] - expected_length) <= max(5, expected_length * 0.2):
                length_consistent_clips += 1

        if length_consistent_clips >= 2:
            insertion_evidence_score += 2
            evidence_details.append(f"Soft clip length consistent with insertion sequence ({length_consistent_clips})")

    # Final decision
    if insertion_evidence_score >= 6:
        conclusion = "Highly likely insertion variant"
        confidence_level = "High"
        is_insertion = True
    elif insertion_evidence_score >= 4:
        conclusion = "Likely insertion variant"
        confidence_level = "Medium"
        is_insertion = True
    elif insertion_evidence_score >= 2:
        conclusion = "Possible insertion variant, further validation needed"
        confidence_level = "Low"
        is_insertion = True
    else:
        conclusion = "Insufficient evidence for insertion"
        confidence_level = "Unsupported"
        is_insertion = False

    return {
        'is_insertion': is_insertion,
        'conclusion': conclusion,
        'confidence_level': confidence_level,
        'evidence_score': insertion_evidence_score,
        'split_reads_count': split_reads,
        'left_clips_count': len(left_clips),
        'right_clips_count': len(right_clips),
        'evidence_details': evidence_details
    }


def contains_insertion_kmers(clip_seq, insertion_kmers, k=21):

    if len(clip_seq) < k or not insertion_kmers:
        return False


    clip_kmers = {clip_seq[i:i + k] for i in range(len(clip_seq) - k + 1)
                  if 'N' not in clip_seq[i:i + k]}

    return len(clip_kmers & insertion_kmers) > 0


def evaluate_germline_evidence_ins(bam_normal, sv_candidate):

    chrom = sv_candidate['start_chr']
    pos = sv_candidate['start_loc']
    inserted_seq = sv_candidate.get('inserted_sequences', 'N/A')

    # 使用插入检测函数
    result = detect_insertion_comprehensive(bam_normal, chrom, pos, inserted_seq)

    if result is None:
        return False, 0, ["Failed"], 0

    evidence_score = result['evidence_score']
    evidence_details = result['evidence_details']
    is_germline = result['is_insertion']
    num_supporting_reads = result['split_reads_count']

    return is_germline, evidence_score, evidence_details, num_supporting_reads



def process_cigar_numba1(cigar_ops, ref_start, min_size, query_sequence_length):

    del_intervals = []
    ins_intervals = []

    reference_pos = ref_start
    query_pos = 0

    for i in range(len(cigar_ops)):
        operation = cigar_ops[i, 1]
        length = cigar_ops[i, 0]

        if operation in (0, 7, 8):  # Match or alignment match
            reference_pos += length
            query_pos += length
        elif operation == 2:  # Deletion
            if length >= min_size:
                del_intervals.append([reference_pos, reference_pos + length, length])
            reference_pos += length
        elif operation == 1:  # Insertion
            if length >= min_size:
                ins_intervals.append([reference_pos, reference_pos, length, query_pos, query_pos + length])
            query_pos += length

    return del_intervals, ins_intervals
def filter_germline_sv_candidates(
    passed_breakpoints,
    bam_file_normal,
    ref_path,
    global_cov
):
    bam_normal = pysam.AlignmentFile(bam_file_normal, "rb")
    ref_file = pysam.FastaFile(ref_path)

    count_yes = 0   # somatic
    count_no = 0    # germline
    countt = 0
    somatic_candidates = []
    global_cov = global_cov * 0.8  # Slightly lower threshold to avoid over-filtering

    for i, candi_filter_short in enumerate(passed_breakpoints):
        chrom_start = candi_filter_short['start_chr']
        chrom_end = candi_filter_short['end_chr']
        start_loc = candi_filter_short['start_loc']
        end_loc = candi_filter_short['end_loc']
        notation = candi_filter_short['breakpoint_notation']

        # Invalid coordinates
        if start_loc < 0 or end_loc < 0:
            print("Invalid coordinates detected, skipping...")
            somatic_candidates.append(candi_filter_short)
            continue

        # Skip large SVs
        if chrom_start == chrom_end:
            length_candi = abs(end_loc - start_loc)
            if length_candi > 50000:
                print("Large SV (>50kb), skipping detailed analysis")
                somatic_candidates.append(candi_filter_short)
                continue

        print(f"\nProcessing candidate {i}: {notation} {chrom_start}:{start_loc}-{chrom_end}:{end_loc}")
        countt += 1

        try:
            is_germline = False
            total_score = 0
            evidence_details = []
            num_kmer_reads = 0

            if chrom_start == chrom_end and notation == '+-':
                print("SV type: DELETION")
                start_loc, end_loc = validate_coordinates(start_loc, end_loc)

                is_kmer_germline, kmer_score, kmer_details, num_kmer_reads = evaluate_germline_evidence(
                    bam_normal, ref_file, candi_filter_short, i, margin=500)

                test_depth = detect_deletion_comprehensive(
                    bam_file=bam_file_normal,
                    chrom=chrom_start,
                    start=start_loc,
                    end=end_loc,
                    global_cov=global_cov,
                    extend_bp=2000,
                    bin_size=10,
                    del_threshold=0.5,
                    flanking_ratio=0.7
                )

                is_depth_germline = test_depth['final_assessment']['is_deletion'] if test_depth else False
                depth_score = test_depth['final_assessment']['confidence_score'] if test_depth else 0

                total_score = kmer_score + depth_score
                evidence_details = kmer_details + [test_depth['final_assessment']['conclusion']] if test_depth else kmer_details
                is_germline = is_kmer_germline or is_depth_germline

            elif chrom_start == chrom_end and notation == '-+':
                print("SV type: DUPLICATION")
                start_loc, end_loc = validate_coordinates(start_loc, end_loc)

                is_kmer_germline, kmer_score, kmer_details, num_kmer_reads = evaluate_germline_evidence_dup(
                    bam_normal, ref_file, candi_filter_short, i, margin=500)

                test_depth = detect_duplication_comprehensive(
                    bam_file=bam_file_normal,
                    chrom=chrom_start,
                    start=start_loc,
                    end=end_loc,
                    global_cov=global_cov,
                    extend_bp=2000,
                    bin_size=10,
                    dup_threshold=1.5,
                    flanking_ratio=1.3
                )

                is_depth_germline = test_depth['final_assessment']['is_duplication'] if test_depth else False
                depth_score = test_depth['final_assessment']['confidence_score'] if test_depth else 0

                total_score = kmer_score + depth_score
                evidence_details = kmer_details + [test_depth['final_assessment']['conclusion']] if test_depth else kmer_details
                is_germline = is_kmer_germline or is_depth_germline

            elif chrom_start == chrom_end and notation in ['++', '--']:
                print("SV type: INVERSION")
                start_loc, end_loc = validate_coordinates(start_loc, end_loc)

                result = detect_inversion_validation(bam_file_normal, chrom_start, start_loc, end_loc)
                is_germline = result['is_inversion']
                total_score = 3 if is_germline else 0
                evidence_details = [result['conclusion']]
                num_kmer_reads = 0

            elif chrom_start != chrom_end:
                print("SV type: TRANSLOCATION")

                result = detect_bnd_by_discordant_pairs(bam_file_normal, chrom_start, start_loc, chrom_end, end_loc)
                is_germline = result['is_bnd']
                total_score = 3 if is_germline else 0
                evidence_details = [result['conclusion']]
                num_kmer_reads = 0

            else:
                print(f"SV type: INSERTION or COMPLEX: {notation}")
                inserted_seq = candi_filter_short.get('inserted_sequences', 'N/A')

                if inserted_seq != 'N/A' and len(inserted_seq) > 5:
                    print("SV type: INSERTION")
                    # print(f"Inserted sequence length: {len(inserted_seq)}bp")

                    is_kmer_germline, kmer_score, kmer_details, num_kmer_reads = evaluate_germline_evidence_ins(
                        bam_file_normal, candi_filter_short)

                    total_score = kmer_score
                    evidence_details = kmer_details
                    is_germline = is_kmer_germline

                else:
                    print("Unknown or complex SV type, conservatively considered somatic")
                    is_germline = False
                    total_score = 0
                    evidence_details = ["Unknown SV type, retained as somatic"]
                    num_kmer_reads = 0

            # Output analysis results
            print(f"  K-mer reads: {num_kmer_reads}")
            print(f"  Evidence score: {total_score}")
            for detail in evidence_details:
                print(f"  Evidence: {detail}")

            if is_germline:
                print(f"GERMLINE DETECTED - Total score: {total_score}")
                print("This SV will be filtered out as germline")
                print('-----------')
                count_no += 1
            else:
                print(f"NO GERMLINE EVIDENCE - Total score: {total_score}")
                print("This SV remains as somatic candidate")
                print('-----------')
                count_yes += 1
                somatic_candidates.append(candi_filter_short)

        except Exception as e:
            print(f"Error processing candidate {i}: {e}")
            print("Adding to somatic candidates to be safe...")
            somatic_candidates.append(candi_filter_short)
            continue

    print(f"\nFiltering summary:")
    print(f"Original candidates: {len(passed_breakpoints)}")
    print(f"Remaining somatic candidates: {len(somatic_candidates)}")
    print(f"Filtered as germline: {len(passed_breakpoints) - len(somatic_candidates)}")

    return somatic_candidates