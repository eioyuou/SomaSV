import pysam

from somasv.utils import validate_coordinates, ensure_int
from somasv.sv_detection import (
    evaluate_germline_evidence, detect_deletion_comprehensive,
    evaluate_germline_evidence_dup, detect_duplication_comprehensive,
    detect_inversion_validation, detect_bnd_by_discordant_pairs,
    evaluate_germline_evidence_ins
)


def filter_germline_sv_candidates(passed_breakpoints, bam_file_normal, ref_path, global_cov):
    bam_normal = pysam.AlignmentFile(bam_file_normal, "rb")
    ref_file = pysam.FastaFile(ref_path)

    count_yes = 0
    count_no = 0
    countt = 0
    somatic_candidates = []
    global_cov = global_cov * 0.8

    for i, candi_filter_short in enumerate(passed_breakpoints):
        chrom_start = candi_filter_short['start_chr']
        chrom_end = candi_filter_short['end_chr']
        start_loc = candi_filter_short['start_loc']
        end_loc = candi_filter_short['end_loc']
        notation = candi_filter_short['breakpoint_notation']

        if start_loc < 0 or end_loc < 0:
            print("Invalid coordinates detected, skipping...")
            somatic_candidates.append(candi_filter_short)
            continue

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
                    bam_file=bam_file_normal, chrom=chrom_start,
                    start=start_loc, end=end_loc, global_cov=global_cov,
                    extend_bp=2000, bin_size=10, del_threshold=0.5, flanking_ratio=0.7
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
                    bam_file=bam_file_normal, chrom=chrom_start,
                    start=start_loc, end=end_loc, global_cov=global_cov,
                    extend_bp=2000, bin_size=10, dup_threshold=1.5, flanking_ratio=1.3
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
                is_germline = result['is_real_inversion']
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