#!/usr/bin/env python3
"""
SomaSV: Long-read Somatic Structural Variant Hunter
Version: 0.0.1
Author: Rt G
"""
import argparse
import os
import sys
import time
import traceback
import subprocess
import warnings
from collections import defaultdict
from somasv import print_banner
warnings.filterwarnings('ignore')

from somasv.utils import get_chrom_lengths
from somasv.multiprocess_runner import run_multiprocess_for_all_chromosomes
from somasv.clustering import parallel_group_related_breakpoints
from somasv.calling import parallel_call_breakpoints
from somasv.depth import compute_depth
from somasv.vcf_io import generate_files, write_somatic_vcf
from somasv.filtering import filter_somatic_variants_improved
from somasv.common_sv import get_common_sv_candidates, identify_somatic_sv
from somasv.germline_filter import filter_germline_sv_candidates


def parse_arguments():
    program_description = """
    SomaSV v0.0.1
    A tool for detecting somatic structural variants using long-read sequencing data.

    Supports two modes:
    1. Long-read only mode: Uses tumor and normal long-read data
    2. Hybrid mode: Uses long-read tumor + normal, plus short-read normal for validation

    For more information, visit: https://github.com/eioyuou/SomaSV
    """
    parser = argparse.ArgumentParser(
        prog='SomaSV',
        description=program_description,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--version', action='version', version='%(prog)s 0.0.1')

    parser.add_argument('--tumor-bam', required=True, help='Path to tumor long-read BAM file')
    parser.add_argument('--normal-bam', required=True, help='Path to normal long-read BAM file')
    parser.add_argument('--reference', required=True, help='Path to reference genome FASTA file')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--sample-platform', required=True, choices=['ONT', 'HIFI', 'ont', 'hifi'],
                        help='Sequencing platform used for long-read samples')

    parser.add_argument('--mode', default='long-read-only', choices=['long-read-only', 'hybrid'],
                        help='Analysis mode (default: long-read-only)')

    parser.add_argument('--tumor-coverage', type=float, help='Tumor sample coverage')
    parser.add_argument('--normal-coverage', type=float, help='Normal sample coverage')

    parser.add_argument('--short-read-normal-bam', help='Path to short-read normal BAM file (hybrid mode)')
    parser.add_argument('--short-read-coverage', type=float, help='Short-read normal sample coverage')
    parser.add_argument('--pon-vcf', help='Path to Panel of Normals (PoN) VCF file for germline filtering (e.g., gnomAD SV sites)')
    parser.add_argument('--short-read-insert-size', type=int, default=250, help='Expected insert size (default: 250)')

    parser.add_argument('--min-length', type=int, default=50, help='Minimum SV length (default: 50)')
    parser.add_argument('--min-mapq', type=int, default=5, help='Minimum mapping quality (default: 5)')
    parser.add_argument('--coverage-binsize', type=int, default=10, help='Coverage bin size (default: 10)')
    parser.add_argument('--internal-num', type=int, default=10000000, help='Internal chunk size (default: 10000000)')

    parser.add_argument('--extension', type=int, default=200, help='Extension size for clustering (default: 200)')
    parser.add_argument('--insertion-additional', type=int, default=250, help='Additional size for insertion clustering (default: 250)')
    parser.add_argument('--end-extension', type=int, default=50, help='Extension size for end clustering (default: 50)')

    parser.add_argument('--min-support', type=int, default=3, help='Minimum support reads (default: 3)')
    parser.add_argument('--bed-file', default=None, help='BED file for additional filtering')

    parser.add_argument('--enable-assembly', action='store_true', help='Enable local assembly (hybrid mode)')
    parser.add_argument('--minimap2-path', default='minimap2', help='Path to minimap2')
    parser.add_argument('--samtools-path', default='samtools', help='Path to samtools')

    parser.add_argument('--processes', type=int, default=20, help='Number of processes (default: 20)')

    parser.add_argument('--debug', action='store_true', help='Enable debug output')
    parser.add_argument('--keep-intermediate', action='store_true', help='Keep intermediate files')

    return parser.parse_args()


def validate_inputs(args):
    """Validate input parameters and files"""
    required_files = [args.tumor_bam, args.normal_bam, args.reference]
    for file_path in required_files:
        if not os.path.exists(file_path):
            sys.exit(f"Error: Required file does not exist: {file_path}")

    for bam_file in [args.tumor_bam, args.normal_bam]:
        if not os.path.exists(bam_file + '.bai'):
            sys.exit(f"Error: BAM index not found for {bam_file}")

    if not os.path.exists(args.reference + '.fai'):
        sys.exit(f"Error: Reference index not found for {args.reference}")

    if args.mode == 'hybrid':
        if not args.short_read_normal_bam:
            sys.exit("Error: --short-read-normal-bam is required for hybrid mode")
        if not os.path.exists(args.short_read_normal_bam):
            sys.exit(f"Error: Short-read normal BAM file does not exist: {args.short_read_normal_bam}")
        if not os.path.exists(args.short_read_normal_bam + '.bai'):
            sys.exit(f"Error: Short-read normal BAM index not found: {args.short_read_normal_bam}.bai")
        if args.pon_vcf:
            if not os.path.exists(args.pon_vcf):
                sys.exit(f"Error: PoN VCF file does not exist: {args.pon_vcf}")
        else:
            print("Warning: No PoN VCF provided. Germline filtering may be less effective.")
        if args.enable_assembly:
            try:
                subprocess.run([args.minimap2_path, '--version'], capture_output=True, check=True)
            except (subprocess.CalledProcessError, FileNotFoundError):
                sys.exit(f"Error: minimap2 not found or not executable: {args.minimap2_path}")
            try:
                subprocess.run([args.samtools_path, '--version'], capture_output=True, check=True)
            except (subprocess.CalledProcessError, FileNotFoundError):
                sys.exit(f"Error: samtools not found or not executable: {args.samtools_path}")

    if args.bed_file and not os.path.exists(args.bed_file):
        sys.exit(f"Error: BED file not found: {args.bed_file}")

    args.sample_platform = args.sample_platform.upper()

    if args.min_length < 1:
        sys.exit("Error: --min-length must be >= 1")
    if args.min_support < 1:
        sys.exit("Error: --min-support must be >= 1")
    if args.processes < 1:
        sys.exit("Error: --processes must be >= 1")

    try:
        os.makedirs(args.output_dir, exist_ok=True)
    except Exception as e:
        sys.exit(f"Error creating output directory: {e}")

    print(f"Validation passed. Running in {args.mode} mode.")
    if args.mode == 'hybrid':
        print(f"Assembly enabled: {args.enable_assembly}")


def main():

    print_banner()
    args = parse_arguments()
    validate_inputs(args)

    contig_order = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']

    chrom_lengths = get_chrom_lengths(args.reference)
    first_24_chrom_lengths = {chrom: chrom_lengths[chrom] for chrom in contig_order if chrom in chrom_lengths}

    print("Starting SomaSV pipeline...")
    print(f"Mode: {args.mode}")
    print(f"Platform: {args.sample_platform}")
    time_start = time.time()

    try:
        # Step 1: Initial SV Detection
        print("Step 1: Running initial SV detection...")
        merged_result, merged_coverage_arrays, global_stats = run_multiprocess_for_all_chromosomes(
            args.tumor_bam, args.normal_bam,
            args.min_length, args.min_mapq,
            contig_order, first_24_chrom_lengths,
            args.internal_num, args.coverage_binsize,
            args.sample_platform, args.processes
        )

        # Step 2: Clustering
        print("Step 2: Clustering breakpoints...")
        clustered_result = parallel_group_related_breakpoints(
            merged_result, args.tumor_bam, args.normal_bam,
            extension=args.extension, insertion_additional=args.insertion_additional,
            num_processes=args.processes
        )

        # Step 3: Breakpoint Calling
        print("Step 3: Calling breakpoints...")
        final_breakpoints = parallel_call_breakpoints(
            clustered_result, args.tumor_bam, args.normal_bam,
            end_extension=args.end_extension, min_length=args.min_length,
            min_support=args.min_support, num_processes=args.processes
        )

        # Step 4: Depth Annotation
        print("Step 4: Annotating breakpoints with depth information...")
        final_annotated_breakpoints = []
        for contig, called_breakpoints in final_breakpoints.items():
            print(f'Annotating breakpoints for contig {contig}...')
            annotated_breakpoints = compute_depth(
                called_breakpoints, merged_coverage_arrays, args.coverage_binsize,
            )
            final_annotated_breakpoints.extend(annotated_breakpoints)

        # Step 5: Generate Initial Output
        print("Step 5: Generating initial output files...")
        generate_files(final_annotated_breakpoints, args.reference, args.output_dir)

        # Step 6: Somatic Filtering
        print("Step 6: Filtering somatic variants...")
        output_prefix = os.path.join(args.output_dir, "filtered")
        final_lengths, passed_breakpoints = filter_somatic_variants_improved(
            final_annotated_breakpoints, output_prefix, args.sample_platform,
            args.reference, global_stats=global_stats, bed_file=args.bed_file,
        )

        print(f"Initial filtering completed: {len(passed_breakpoints)} variants passed")

        # Mode-specific processing
        if args.mode == 'long-read-only':
            final_somatic_candidates = passed_breakpoints
            print("Long-read only mode: Analysis completed")

        elif args.mode == 'hybrid':
            print("Step 7: Hybrid mode processing...")

            # Step 7a: Common SV Filtering
            if args.pon_vcf:
                print("Step 7a: Filtering against Panel of Normals (PoN)...")
                sv_by_chrom_long = defaultdict(list)
                for bp in passed_breakpoints:
                    chrom = bp['start_chr']
                    sv_by_chrom_long[chrom].append(bp)

                common_sv_candidates = get_common_sv_candidates(args.pon_vcf)
                after_common_filtering = identify_somatic_sv(sv_by_chrom_long, common_sv_candidates)
                print(f"After common SV filtering: {len(after_common_filtering)} variants remain")
            else:
                after_common_filtering = passed_breakpoints
                print("No PoN VCF provided, skipping PoN filtering")

            # Step 7b: Germline Evidence Filtering (Short-read)
            print("Step 7b: Germline evidence filtering with short reads...")

            if args.short_read_coverage is None:
                print("Calculating short-read coverage...")
                estimated_coverage = 30
                print(f"Using estimated coverage: {estimated_coverage}x")
            else:
                estimated_coverage = args.short_read_coverage

            after_shortread_filtering = filter_germline_sv_candidates(
                passed_breakpoints=after_common_filtering,
                bam_file_normal=args.short_read_normal_bam,
                ref_path=args.reference,
                global_cov=estimated_coverage
            )
            print(f"After short-read filtering: {len(after_shortread_filtering)} variants remain")

            # Output intermediate results
            print("Generating intermediate output (after short-read filtering)...")
            intermediate_vcf_path = os.path.join(args.output_dir, "shortread_filtered_somatic_variants.vcf")
            write_somatic_vcf(after_shortread_filtering, intermediate_vcf_path, args.reference)
            print(f"Intermediate results saved to: {intermediate_vcf_path}")

            final_somatic_candidates = after_shortread_filtering

        # Step 8: Generate Final Output
        print("Step 8: Generating final output...")
        final_vcf_path = os.path.join(args.output_dir, "final_somatic_variants.vcf")
        write_somatic_vcf(final_somatic_candidates, final_vcf_path, args.reference)

        # Summary
        print("\n" + "=" * 60)
        print("SomaSV Analysis Summary")
        print("=" * 60)
        print(f"Mode: {args.mode}")
        print(f"Platform: {args.sample_platform}")
        print(f"Total runtime: {time.time() - time_start:.2f} seconds")
        print(f"Initial breakpoints detected: {len(final_annotated_breakpoints)}")
        print(f"After initial filtering: {len(passed_breakpoints)}")
        if args.mode == 'hybrid':
            if args.pon_vcf:
                print(f"After PoN filtering: {len(after_common_filtering)}")
            print(f"After short-read filtering: {len(after_shortread_filtering)}")
            print(f"  -> Intermediate VCF: {intermediate_vcf_path}")
        print(f"Final somatic variants: {len(final_somatic_candidates)}")
        print(f"Final VCF: {final_vcf_path}")
        print("=" * 60)
        print("Analysis completed successfully!")

    except Exception as e:
        print(f"Error during analysis: {str(e)}")
        if args.debug:
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()