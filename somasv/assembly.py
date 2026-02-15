import os
import re
import pysam
import multiprocessing
import tempfile
import shutil

from mantaAss_module import mantaAss
from somasv.short_read_analysis import (
    is_read_filtered_short, extract_split_reads_short,
    get_sv_candidates_from_read_indels_short, extract_discordant_pairs_short,
    extract_semi_aligned_short, collect_evidence_reads_small_short,
    collect_evidence_reads_large_short, collect_evidence_reads_BND_short
)


def save_contig_to_fastq(contig_seq, output_fastq, contig_id="assembled_contig"):
    with open(output_fastq, "a") as fq:
        fq.write(f"@{contig_id}\n")
        fq.write(f"{contig_seq}\n")
        fq.write("+\n")
        fq.write("I" * len(contig_seq) + "\n")


def has_assemble_in_region(candidates, count, save_path, sv_type):
    if not candidates:
        print(f"No reads for assembly in {sv_type}")
        return False, "No reads for assembly"

    try:
        min_k, max_k = 35, 75
        contig_num = 10
        assembled_read_lis = mantaAss(candidates, min_k, max_k, contig_num)

        contig_count = 0
        for i, contig_seq in enumerate(assembled_read_lis[1]):
            if not contig_seq:
                continue
            contig_name = f'contig_{count}_{sv_type}_{i}'

            if isinstance(contig_seq, bytes):
                contig_seq_str = contig_seq.decode()
            else:
                contig_seq_str = str(contig_seq)

            save_contig_to_fastq(contig_seq_str, save_path, contig_id=contig_name)
            contig_count += 1

        print(f"Assembled {contig_count} contigs for {sv_type}")
        return True, f"Assembled {contig_count} contigs"

    except Exception as e:
        print(f"Assembly failed for {sv_type}: {e}")
        return False, f"Assembly failed: {e}"


def process_single_sv_task(args):
    (i, candi_filter_short, bam_file_path, ref_fasta_path,
     sample_name, min_mapq, base_isize, length_candi_threshold,
     min_indel, max_isize, chrom_list, ref_chr_length,
     ref_chr_id_to_name, output_dir) = args

    bam_normal = pysam.AlignmentFile(bam_file_path, "rb")
    ref_file = pysam.FastaFile(ref_fasta_path)

    chrom_start = candi_filter_short['start_chr']
    chrom_end = candi_filter_short['end_chr']
    start_loc = candi_filter_short['start_loc']
    end_loc = candi_filter_short['end_loc']
    notation = candi_filter_short['breakpoint_notation']
    length_candi = abs(end_loc - start_loc)

    print(f"[#{i}] Processing SV: {notation} {chrom_start}:{start_loc}-{chrom_end}:{end_loc}")

    local_fastq_path = os.path.join(output_dir, f"contig_{i}.fastq")

    try:
        if chrom_start != chrom_end:
            candi_start, candi_end = collect_evidence_reads_BND_short(
                bam_normal, chrom_start, chrom_end, start_loc, end_loc,
                min_mapq, ref_chr_length, sample_name, chrom_list,
                ref_chr_id_to_name, ref_fasta_path,
                min_isize_large=base_isize + 500,
                max_isize_large=max_isize,
                breakpoint_margin=1000
            )
            has_assemble_in_region(candi_start, i, local_fastq_path, 'BND_start')
            has_assemble_in_region(candi_end, i, local_fastq_path, 'BND_end')
        else:
            if notation == '+-':
                if length_candi < length_candi_threshold:
                    candi_read = collect_evidence_reads_small_short(
                        bam_normal, chrom_start, start_loc, end_loc,
                        min_mapq, ref_chr_length, sample_name, chrom_list,
                        ref_chr_id_to_name, ref_fasta_path,
                        min_indel, max_isize,
                        breakpoint_margin=500 if length_candi < 500 else (length_candi + 500)
                    )
                    has_assemble_in_region(candi_read, i, local_fastq_path, 'small_DEL')
                else:
                    candi_start, candi_end = collect_evidence_reads_large_short(
                        bam_normal, chrom_start, start_loc, end_loc,
                        min_mapq, ref_chr_length, sample_name, chrom_list,
                        ref_chr_id_to_name, ref_fasta_path,
                        min_isize_large=base_isize + 500,
                        max_isize_large=max_isize,
                        breakpoint_margin=1000
                    )
                    has_assemble_in_region(candi_start, i, local_fastq_path, 'DEL_start')
                    has_assemble_in_region(candi_end, i, local_fastq_path, 'DEL_end')
            elif notation == '-+':
                if length_candi < length_candi_threshold:
                    candi_read = collect_evidence_reads_small_short(
                        bam_normal, chrom_start, start_loc, end_loc,
                        min_mapq, ref_chr_length, sample_name, chrom_list,
                        ref_chr_id_to_name, ref_fasta_path,
                        min_indel, max_isize,
                        breakpoint_margin=500 if length_candi < 500 else (length_candi + 500)
                    )
                    has_assemble_in_region(candi_read, i, local_fastq_path, 'small_DUP')
                else:
                    candi_start, candi_end = collect_evidence_reads_large_short(
                        bam_normal, chrom_start, start_loc, end_loc,
                        min_mapq, ref_chr_length, sample_name, chrom_list,
                        ref_chr_id_to_name, ref_fasta_path,
                        min_isize_large=base_isize + 500,
                        max_isize_large=max_isize,
                        breakpoint_margin=1000
                    )
                    has_assemble_in_region(candi_start, i, local_fastq_path, 'DUP_start')
                    has_assemble_in_region(candi_end, i, local_fastq_path, 'DUP_end')
            elif notation in ['++', '--']:
                if length_candi < length_candi_threshold:
                    candi_read = collect_evidence_reads_small_short(
                        bam_normal, chrom_start, start_loc, end_loc,
                        min_mapq, ref_chr_length, sample_name, chrom_list,
                        ref_chr_id_to_name, ref_fasta_path,
                        min_indel, max_isize,
                        breakpoint_margin=500 if length_candi < 500 else (length_candi + 500)
                    )
                    has_assemble_in_region(candi_read, i, local_fastq_path, 'small_INV')
                else:
                    candi_start, candi_end = collect_evidence_reads_large_short(
                        bam_normal, chrom_start, start_loc, end_loc,
                        min_mapq, ref_chr_length, sample_name, chrom_list,
                        ref_chr_id_to_name, ref_fasta_path,
                        min_isize_large=base_isize + 500,
                        max_isize_large=max_isize,
                        breakpoint_margin=1000
                    )
                    has_assemble_in_region(candi_start, i, local_fastq_path, 'INV_start')
                    has_assemble_in_region(candi_end, i, local_fastq_path, 'INV_end')

    except Exception as e:
        print(f"[#{i}] Error while processing SV: {e}")
    finally:
        bam_normal.close()
        ref_file.close()

    return local_fastq_path


def process_sv_list_and_assemble_contigs_short_mp(
    somatic_sv_list, bam_file_path, ref_fasta_path, output_fastq_path,
    sample_name='normal', min_mapq=0, base_isize=250,
    length_candi_threshold=10000, min_indel=40, max_isize=1000000, num_workers=20
):
    chrom_list = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']

    ref_file = pysam.FastaFile(ref_fasta_path)
    ref_chr_length = {
        ref: ref_file.get_reference_length(ref)
        for ref in chrom_list if ref in ref_file.references
    }

    bam_normal = pysam.AlignmentFile(bam_file_path, "rb")
    ref_chr_id_to_name = {
        i: name for i, name in enumerate(bam_normal.references)
    }

    temp_dir = tempfile.mkdtemp(prefix='sv_assemble_')

    param_list = []
    for i, sv in enumerate(somatic_sv_list):
        param_list.append((
            i, sv, bam_file_path, ref_fasta_path,
            sample_name, min_mapq, base_isize, length_candi_threshold,
            min_indel, max_isize,
            chrom_list, ref_chr_length, ref_chr_id_to_name,
            temp_dir
        ))

    if num_workers is None:
        num_workers = max(1, multiprocessing.cpu_count() - 1)
    with multiprocessing.Pool(processes=num_workers) as pool:
        result_fastq_paths = pool.map(process_single_sv_task, param_list)

    with open(output_fastq_path, 'w') as out_f:
        for fq in result_fastq_paths:
            if os.path.exists(fq):
                with open(fq) as f:
                    out_f.write(f.read())

    shutil.rmtree(temp_dir)
    print(f"\nMulti-process SV assembly completed. Output: {output_fastq_path}")


def align_contigs_to_reference(
    contig_fastq, ref_path, output_prefix,
    minimap2_path="minimap2", samtools_path="samtools", preset="asm5"
):
    import subprocess
    sam_path = f"{output_prefix}.sam"
    bam_path = f"{output_prefix}.bam"
    sorted_bam_path = f"{output_prefix}_sorted.bam"

    print(f"Step 1: Running minimap2...")
    minimap2_cmd = [minimap2_path, "-ax", preset, ref_path, contig_fastq]
    with open(sam_path, "w") as sam_file:
        subprocess.run(minimap2_cmd, stdout=sam_file, check=True)

    print(f"Step 2: Converting SAM to BAM...")
    subprocess.run([samtools_path, "view", "-bS", sam_path, "-o", bam_path], check=True)

    print(f"Step 3: Sorting BAM file...")
    subprocess.run([samtools_path, "sort", "-o", sorted_bam_path, bam_path], check=True)

    print(f"Step 4: Indexing sorted BAM file...")
    subprocess.run([samtools_path, "index", sorted_bam_path], check=True)

    print(f"Alignment completed. Sorted BAM: {sorted_bam_path}")
    return sorted_bam_path


def filter_unmatched_breakpoints(somatic_sv_list, potential_breakpoints_contig, distance_threshold=10000):
    unmatched_breakpoints = []

    for i, match_bp in enumerate(somatic_sv_list):
        matched = False
        for bp in potential_breakpoints_contig:
            match = re.search(r"contig_(\d+)_", bp['read_name'])
            if not match:
                continue
            number = int(match.group(1))
            if number == i:
                if (
                    bp['start_chr'] == match_bp['start_chr']
                    and bp['end_chr'] == match_bp['end_chr']
                    and abs(bp['start_loc'] - match_bp['start_loc']) <= distance_threshold
                    and abs(bp['end_loc'] - match_bp['end_loc']) <= distance_threshold
                ):
                    matched = True
                    break
        if not matched:
            unmatched_breakpoints.append(match_bp)

    return unmatched_breakpoints