import os
import csv
import re
import pysam
from collections import defaultdict
from intervaltree import Interval, IntervalTree


# ==================== BED Support ====================

def get_support_read(chrom, pos, bed_trees):
    if not bed_trees or chrom not in bed_trees:
        return None
    overlaps = list(bed_trees[chrom].overlap(pos - 1, pos + 1))
    if not overlaps:
        return None
    return overlaps[0].data


def load_bed_regions(bed_file):
    if not bed_file or not os.path.exists(bed_file):
        return None

    bed_trees = defaultdict(IntervalTree)

    try:
        with open(bed_file) as f:
            header_lines = 0
            for line in f:
                if line.startswith('#'):
                    header_lines += 1
                    continue

                if header_lines < 2:
                    print("Warning: Bed file might be missing header lines")

                fields = line.strip().split('\t')
                if len(fields) < 11:
                    continue

                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2]) + 1
                support_info = {
                    'cluster_size': int(fields[3]),
                    'mean_pos': float(fields[4]),
                    'std_dev': float(fields[5]),
                    'mean_quality': float(fields[6]),
                    'support_reads_tumor': int(fields[7]),
                    'support_reads_normal': int(fields[8]),
                    'strand_ratio_tumor': float(fields[9]),
                    'strand_ratio_normal': float(fields[10])
                }

                bed_trees[chrom][start:end] = support_info
    except Exception as e:
        print(f"Error loading BED file: {e}")
        return None

    return bed_trees


# ==================== VCF Reading ====================

def get_info_fields(vcf_file):
    fields = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('##INFO='):
                info_str = line.strip()[8:-1]
                parts = info_str.split(',')
                for part in parts:
                    if part.startswith('ID='):
                        fields.append(part[3:])
                        break
            elif line.startswith('#CHROM'):
                break
    return fields


def parse_info_field(info_str):
    info_dict = {}
    if not info_str or info_str == '.':
        return info_dict
    items = info_str.split(';')
    for item in items:
        if '=' in item:
            key, value = item.split('=', 1)
            try:
                if ',' in value:
                    parts = value.split(',')
                    try:
                        value = [float(part) for part in parts]
                    except ValueError:
                        pass
                else:
                    try:
                        value = float(value)
                    except ValueError:
                        pass
            except:
                pass
            info_dict[key] = value
        else:
            info_dict[item] = True
    return info_dict


def read_vcf(vcf_file):
    header = []
    data = []
    samples = []
    print(f'Reading {vcf_file}')

    try:
        with open(vcf_file, 'r') as f:
            info_fields = []
            for line in f:
                if line.startswith('##INFO='):
                    info_str = line.strip()[8:-1]
                    parts = info_str.split(',')
                    for part in parts:
                        if part.startswith('ID='):
                            info_fields.append(part[3:])
                            break
                elif line.startswith('#CHROM'):
                    columns = line.strip().split('\t')
                    samples = columns[9:]
                    header = ['SAMPLE', 'ID', 'CHROM', 'POS', 'CHR2', 'END', 'SVTYPE']
                    header.extend(info_fields)
                    break

            if len(samples) != 1:
                raise ValueError(
                    f'One sample per VCF permitted: {len(samples)} samples found in {vcf_file} ({",".join(samples)})')

            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 8:
                    continue

                chrom = fields[0]
                pos = int(fields[1])
                var_id = fields[2]
                ref = fields[3]
                alt = fields[4]
                qual = fields[5]
                filter_str = fields[6]
                info_str = fields[7]

                info_dict = parse_info_field(info_str)

                row = {
                    'SAMPLE': samples[0] if samples else 'UNKNOWN',
                    'ID': var_id,
                    'CHROM': chrom,
                    'POS': pos,
                    'REF': ref,
                    'ALT': alt,
                    'QUAL': qual,
                    'FILTER': filter_str,
                    'INFO': info_str,
                    'CHR2': info_dict.get('CHR2', chrom),
                    'END': info_dict.get('END', pos),
                    'SVLEN': info_dict.get('SVLEN', 0),
                    'SVTYPE': info_dict.get('SVTYPE', '')
                }

                for field in info_fields:
                    if field in info_dict:
                        row[field] = info_dict[field]
                    else:
                        row[field] = None

                data.append(row)

    except Exception as e:
        print(f"Error reading VCF file: {e}")
        return header, []

    print(f'Done reading {vcf_file}, found {len(data)} variants')
    return header, data


def parse_info(info_str):
    info_dict = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[item] = True
    return info_dict


def format_data(data_matrix):
    updated_data_matrix = []

    for row_data in data_matrix:
        row_data['POS'] = int(row_data['POS']) if row_data['POS'] else 0
        row_data['END'] = int(row_data['END']) if row_data['END'] else 0

        for dp_field in ['TUMOR_DP_BEFORE', 'TUMOR_DP_AT', 'TUMOR_DP_AFTER', 'NORMAL_DP_BEFORE', 'NORMAL_DP_AT',
                         'NORMAL_DP_AFTER']:
            if dp_field in row_data:
                if isinstance(row_data[dp_field], (tuple, list)):
                    dp_values = row_data[dp_field]
                else:
                    dp_values = row_data[dp_field].split(',') if row_data[dp_field] else ['0', '0']

                row_data[f'{dp_field}_0'] = dp_values[0]
                row_data[f'{dp_field}_1'] = dp_values[1] if len(dp_values) > 1 else dp_values[0]

        for af_field in ['TUMOR_AF', 'NORMAL_AF']:
            if af_field in row_data:
                if isinstance(row_data[af_field], (tuple, list)):
                    af_values = row_data[af_field]
                else:
                    af_values = row_data[af_field].split(',') if row_data[af_field] else ['0', '0']

                row_data[f'{af_field}_0'] = af_values[0]
                row_data[f'{af_field}_1'] = af_values[1] if len(af_values) > 1 else af_values[0]

        for key in ['TUMOR_DP_BEFORE_1', 'TUMOR_DP_AT_1', 'TUMOR_DP_AFTER_1', 'NORMAL_DP_BEFORE_1', 'NORMAL_DP_AT_1',
                     'NORMAL_DP_AFTER_1', 'TUMOR_AF_1', 'NORMAL_AF_1']:
            if key in row_data and (row_data[key] == '' or row_data[key] is None):
                row_data[key] = row_data[key.replace('_1', '_0')]

        for key, value in row_data.items():
            if value == 'inf' or value == '-inf':
                row_data[key] = '-1'

        if 'ORIGIN_STARTS_STD_DEV' in row_data and 'ORIGIN_EVENT_SIZE_MEAN' in row_data:
            row_data['ORIGIN_STD_MEAN_RATIO'] = float(row_data['ORIGIN_STARTS_STD_DEV']) / (
                    float(row_data['ORIGIN_EVENT_SIZE_MEAN']) + 1.0)
        if 'END_STARTS_STD_DEV' in row_data and 'END_EVENT_SIZE_MEAN' in row_data:
            row_data['END_STD_MEAN_RATIO'] = float(row_data['END_STARTS_STD_DEV']) / (
                    float(row_data['END_EVENT_SIZE_MEAN']) + 1.0)

        row_data['ORIGIN_STD_MEAN_RATIO'] = row_data.get('ORIGIN_STD_MEAN_RATIO', 0)
        row_data['END_STD_MEAN_RATIO'] = row_data.get('END_STD_MEAN_RATIO', 0)

        if 'SVTYPE' in row_data:
            row_data['SVTYPE_NUM'] = {'BND': 0, 'INS': 1, 'SBND': 2}.get(row_data['SVTYPE'], -1)

        updated_data_matrix.append(row_data)

    return updated_data_matrix


# ==================== VCF Writing ====================

def get_stats_str(consensus_breakpoint):
    stats_str = ''
    if 'start_cluster' in consensus_breakpoint and 'stats' in consensus_breakpoint['start_cluster']:
        stats_originating = consensus_breakpoint['start_cluster']['stats']
        for key, value in stats_originating.items():
            stats_str += f'ORIGIN_{key.upper()}={value};'
    if 'end_cluster' in consensus_breakpoint and 'stats' in consensus_breakpoint['end_cluster']:
        stats_end = consensus_breakpoint['end_cluster']['stats']
        for key, value in stats_end.items():
            stats_str += f'END_{key.upper()}={value};'
    return stats_str


def combine_hp_values(values):
    combined_values = []
    current_value = ""
    for v in values:
        if v.isdigit():
            current_value += v
        else:
            if current_value:
                combined_values.append(current_value)
                current_value = ""
            combined_values.append(v)
    if current_value:
        combined_values.append(current_value)
    return [v for v in combined_values if v != ',']


def get_alts(consensus_breakpoint, start_base, end_base):
    if consensus_breakpoint['breakpoint_notation'] == "<INS>":
        return ["<INS>"]
    if consensus_breakpoint['breakpoint_notation'] in ["-", "+"]:
        return [f'.{end_base}']

    alts = ['', '']
    if consensus_breakpoint['breakpoint_notation'].startswith("+"):
        alts[0] += f'{start_base}'
        if consensus_breakpoint['breakpoint_notation'].endswith("+"):
            alts[0] += f']{consensus_breakpoint["end_chr"]}:{consensus_breakpoint["end_loc"]}]'
            alts[1] += f'{end_base}]{consensus_breakpoint["start_chr"]}:{consensus_breakpoint["start_loc"]}]'
        else:
            alts[0] += f'[{consensus_breakpoint["end_chr"]}:{consensus_breakpoint["end_loc"]}['
            alts[1] += f']{consensus_breakpoint["start_chr"]}:{consensus_breakpoint["start_loc"]}]{end_base}'
    else:
        if consensus_breakpoint['breakpoint_notation'].endswith("+"):
            alts[0] += f']{consensus_breakpoint["end_chr"]}:{consensus_breakpoint["end_loc"]}]'
            alts[1] += f'{end_base}[{consensus_breakpoint["start_chr"]}:{consensus_breakpoint["start_loc"]}['
        else:
            alts[0] += f'[{consensus_breakpoint["end_chr"]}:{consensus_breakpoint["end_loc"]}['
            alts[1] += f'[{consensus_breakpoint["start_chr"]}:{consensus_breakpoint["start_loc"]}[{end_base}'
        alts[0] += f'{start_base}'

    return alts


def generate_vcf_header(reference_fasta_path, reference_fai_path, sample_name):
    vcf_header_str = []
    vcf_header_str.extend([
        "##fileformat=VCFv4.2",
        "##source=SomaticHunter0.0.1"
    ])

    assembly_name = os.path.basename(reference_fasta_path)
    with open(reference_fai_path) as f:
        reader = csv.reader(f, delimiter='\t')
        for line in reader:
            contig = line[0]
            length = line[1]
            vcf_header_str.append(f'##contig=<ID={contig},length={length},assembly={assembly_name}>')

    cmd_string = '##SomaticHunter_args="'
    cmd_string += '"'

    vcf_header_str.extend([
        cmd_string,
        f'##reference={{{reference_fasta_path}}}',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        '##INFO=<ID=MATEID,Number=1,Type=String,Description="ID of mate breakends">',
        '##INFO=<ID=NORMAL_READ_SUPPORT,Number=1,Type=Integer,Description="Number of SV supporting normal reads">',
        '##INFO=<ID=TUMOR_READ_SUPPORT,Number=1,Type=Integer,Description="Number of SV supporting tumor reads">',
        '##INFO=<ID=NORMAL_ALN_SUPPORT,Number=1,Type=Integer,Description="Number of SV supporting normal alignments">',
        '##INFO=<ID=TUMOR_ALN_SUPPORT,Number=1,Type=Integer,Description="Number of SV supporting tumor alignments">',
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">',
        '##INFO=<ID=TUMOR_DP_BEFORE,Number=2,Type=Integer,Description="Local tumor depth in bin before the breakpoint(s) of an SV">',
        '##INFO=<ID=TUMOR_DP_AT,Number=2,Type=Integer,Description="Local tumor depth in bin at the breakpoint(s) of an SV">',
        '##INFO=<ID=TUMOR_DP_AFTER,Number=2,Type=Integer,Description="Local tumor depth in bin after the breakpoint(s) of an SV">',
        '##INFO=<ID=NORMAL_DP_BEFORE,Number=2,Type=Integer,Description="Local normal depth in bin before the breakpoint(s) of an SV">',
        '##INFO=<ID=NORMAL_DP_AT,Number=2,Type=Integer,Description="Local normal depth in bin at the breakpoint(s) of an SV">',
        '##INFO=<ID=NORMAL_DP_AFTER,Number=2,Type=Integer,Description="Local normal depth in bin after the breakpoint(s) of an SV">',
        '##INFO=<ID=TUMOR_AF,Number=2,Type=Float,Description="Allele-fraction (AF) of tumor variant-supporting reads to tumor read depth (DP) at breakpoint">',
        '##INFO=<ID=NORMAL_AF,Number=2,Type=Float,Description="Allele-fraction (AF) of normal variant-supporting reads to normal read depth (DP) at breakpoint">',
        '##INFO=<ID=BP_NOTATION,Number=1,Type=String,Description="+- notation format of variant (same for paired breakpoints)">',
        '##INFO=<ID=SOURCE,Number=1,Type=String,Description="Source of evidence for a breakpoint - CIGAR (INS, DEL, SOFTCLIP), SUPPLEMENTARY or mixture">',
        '##INFO=<ID=CLUSTERED_READS_TUMOR,Number=1,Type=Integer,Description="Total number of tumor reads clustered at this location of any SV type">',
        '##INFO=<ID=CLUSTERED_READS_NORMAL,Number=1,Type=Integer,Description="Total number of normal reads clustered at this location of any SV type">',
        '##INFO=<ID=TUMOR_ALT_HP,Number=3,Type=Integer,Description="Counts of SV-supporting reads belonging to each haplotype in the tumor sample (1/2/NA)">',
        '##INFO=<ID=TUMOR_PS,Number=.,Type=String,Description="List of unique phase sets from the tumor supporting reads">',
        '##INFO=<ID=NORMAL_ALT_HP,Number=3,Type=Integer,Description="Counts of reads belonging to each haplotype in the normal sample (1/2/NA)">',
        '##INFO=<ID=NORMAL_PS,Number=.,Type=String,Description="List of unique phase sets from the normal supporting reads">',
        '##INFO=<ID=TUMOR_TOTAL_HP_AT,Number=3,Type=Integer,Description="Counts of all reads at SV location belonging to each haplotype in the tumor sample (1/2/NA)">',
        '##INFO=<ID=NORMAL_TOTAL_HP_AT,Number=3,Type=Integer,Description="Counts of all reads at SV location belonging to each haplotype in the normal sample (1/2/NA)">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">',
        '##INFO=<ID=ORIGIN_STARTS_STD_DEV,Number=1,Type=Float,Description="Originating cluster value for starts_std_dev">',
        '##INFO=<ID=ORIGIN_MAPQ_MEAN,Number=1,Type=Float,Description="Originating cluster value for mapq_mean">',
        '##INFO=<ID=ORIGIN_EVENT_SIZE_STD_DEV,Number=1,Type=Float,Description="Originating cluster value for event_size_std_dev">',
        '##INFO=<ID=ORIGIN_EVENT_SIZE_MEDIAN,Number=1,Type=Float,Description="Originating cluster value for event_size_median">',
        '##INFO=<ID=ORIGIN_EVENT_SIZE_MEAN,Number=1,Type=Float,Description="Originating cluster value for event_size_mean">',
        '##INFO=<ID=END_STARTS_STD_DEV,Number=1,Type=Float,Description="End cluster value for starts_std_dev">',
        '##INFO=<ID=END_MAPQ_MEAN,Number=1,Type=Float,Description="End cluster value for mapq_mean">',
        '##INFO=<ID=END_EVENT_SIZE_STD_DEV,Number=1,Type=Float,Description="End cluster value for event_size_std_dev">',
        '##INFO=<ID=END_EVENT_SIZE_MEDIAN,Number=1,Type=Float,Description="End cluster value for event_size_median">',
        '##INFO=<ID=END_EVENT_SIZE_MEAN,Number=1,Type=Float,Description="End cluster value for event_size_mean">'
    ])

    vcf_header_str.append(
        "#" + "\t".join(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sample_name]))

    return "\n".join(vcf_header_str) + "\n"


def generate_vcf_line(breakpoint_dict, ref_fasta):
    try:
        start_base = ref_fasta.fetch(breakpoint_dict['start_chr'], breakpoint_dict['start_loc'] - 1,
                                     breakpoint_dict['start_loc'])
        start_base = "N" if start_base == "" else start_base
        if start_base.upper() not in ['N', 'A', 'T', 'C', 'G']:
            start_base = "N"
    except ValueError:
        start_base = 'N'

    try:
        end_base = ref_fasta.fetch(breakpoint_dict['end_chr'], breakpoint_dict['end_loc'] - 1,
                                   breakpoint_dict['end_loc'])
        end_base = "N" if end_base == "" else end_base
        if end_base.upper() not in ['N', 'A', 'T', 'C', 'G']:
            end_base = 'N'
    except ValueError:
        end_base = 'N'

    alts = get_alts(breakpoint_dict, start_base, end_base)
    gt_tag = '0/0' if breakpoint_dict['read_support_counts'].get('normal', 0) >= 1 else '0/1'
    info = []

    if breakpoint_dict['breakpoint_notation'] == "<INS>":
        svtype = "INS"
    elif breakpoint_dict['breakpoint_notation'] in ["+", "-"]:
        svtype = "SBND"
    else:
        svtype = "BND"

    info.append(f'SVTYPE={svtype}')
    info.append(f'NORMAL_READ_SUPPORT={breakpoint_dict["read_support_counts"].get("normal", 0)}')
    info.append(f'TUMOR_READ_SUPPORT={breakpoint_dict["read_support_counts"].get("tumor", 0)}')
    info.append(f'NORMAL_ALN_SUPPORT={breakpoint_dict["aln_support_counts"].get("normal", 0)}')
    info.append(f'TUMOR_ALN_SUPPORT={breakpoint_dict["aln_support_counts"].get("tumor", 0)}')

    svlen = breakpoint_dict.get("sv_length", 0)
    info.append(f'SVLEN={svlen}')

    if 'region_total_depths' in breakpoint_dict:
        for sample in ['tumor', 'normal']:
            if sample in breakpoint_dict['region_total_depths']:
                depths = breakpoint_dict['region_total_depths'][sample]
                before_depths = f"{depths[0][0]},{depths[0][1]}"
                at_depths = f"{depths[1][0]},{depths[1][1]}"
                after_depths = f"{depths[2][0]},{depths[2][1]}"
                info.append(f'{sample.upper()}_DP_BEFORE={before_depths}')
                info.append(f'{sample.upper()}_DP_AT={at_depths}')
                info.append(f'{sample.upper()}_DP_AFTER={after_depths}')

    if 'allele_fractions' in breakpoint_dict:
        for sample in ['tumor', 'normal']:
            if sample in breakpoint_dict['allele_fractions']:
                af_values = breakpoint_dict['allele_fractions'][sample]
                af_str = ",".join([str(v if v is not None else 0) for v in af_values])
                info.append(f'{sample.upper()}_AF={af_str}')

    info.append(f'BP_NOTATION={breakpoint_dict.get("breakpoint_notation", "0")}')
    info.append(f'SOURCE={breakpoint_dict.get("source", "0")}')

    tumor_clustered_reads = breakpoint_dict['total_read_counts'].get('tumor', 0)
    normal_clustered_reads = breakpoint_dict['total_read_counts'].get('normal', 0)
    info.append(f'CLUSTERED_READS_TUMOR={tumor_clustered_reads}')
    info.append(f'CLUSTERED_READS_NORMAL={normal_clustered_reads}')

    if 'TUMOR_ALT_HP' in breakpoint_dict:
        info.append(f'TUMOR_ALT_HP={breakpoint_dict["TUMOR_ALT_HP"]}')
    if 'NORMAL_ALT_HP' in breakpoint_dict:
        info.append(f'NORMAL_ALT_HP={breakpoint_dict["NORMAL_ALT_HP"]}')
    if 'TUMOR_PS' in breakpoint_dict:
        info.append(f'TUMOR_PS={breakpoint_dict["TUMOR_PS"]}')
    if 'NORMAL_PS' in breakpoint_dict:
        info.append(f'NORMAL_PS={breakpoint_dict["NORMAL_PS"]}')
    if 'TUMOR_TOTAL_HP_AT' in breakpoint_dict:
        info.append(f'TUMOR_TOTAL_HP_AT={breakpoint_dict["TUMOR_TOTAL_HP_AT"]}')
    if 'NORMAL_TOTAL_HP_AT' in breakpoint_dict:
        info.append(f'NORMAL_TOTAL_HP_AT={breakpoint_dict["NORMAL_TOTAL_HP_AT"]}')

    stats_str = get_stats_str(breakpoint_dict)
    if stats_str:
        info.append(stats_str.rstrip(';'))

    vcf_lines = [[
        breakpoint_dict['start_chr'],
        str(breakpoint_dict['start_loc']),
        f'ID_{breakpoint_dict["count"]}_1',
        start_base,
        alts[0],
        '.',
        'PASS',
        ';'.join([item for item in info if item]),
        'GT',
        gt_tag
    ]]

    if breakpoint_dict['breakpoint_notation'] not in ["<INS>", "+", "-"]:
        vcf_lines.append([
            breakpoint_dict['end_chr'],
            str(breakpoint_dict['end_loc']),
            f'ID_{breakpoint_dict["count"]}_2',
            end_base,
            alts[1],
            '.',
            'PASS',
            ';'.join([item for item in info if item]),
            'GT',
            gt_tag
        ])

    vcf_string = ''
    for line in vcf_lines:
        vcf_string += "\t".join(line) + "\n"

    return vcf_string


def generate_files(annotated_breakpoints, ref_fasta_path, output_dir):
    vcf_file = os.path.join(output_dir, 'all_sv_breakpoints.vcf')
    vcf_file_normal = os.path.join(output_dir, 'all_sv_breakpoints_normal.vcf')

    ref_fasta = pysam.FastaFile(ref_fasta_path)

    vcf_string = generate_vcf_header(ref_fasta_path, ref_fasta_path + '.fai', 'sample_name')
    vcf_string_normal = generate_vcf_header(ref_fasta_path, ref_fasta_path + '.fai', 'sample_name')
    count = 0

    for bp in annotated_breakpoints:
        bp['count'] = count
        vcf_string += generate_vcf_line(bp, ref_fasta)
        count += 1

    with open(vcf_file, 'w') as vcf_output:
        vcf_output.write(vcf_string)

    with open(vcf_file_normal, 'w') as vcf_output_normal:
        vcf_output_normal.write(vcf_string_normal)

    print(f"Files generated:")
    print(f"- VCF: {vcf_file}")
    print(f"- Normal VCF: {vcf_file_normal}")
    print(f"Total breakpoints processed: {count}")


def generate_somatic_vcf_header(ref_fasta_path, sample_name='sample_name'):
    vcf_header = generate_vcf_header(ref_fasta_path, ref_fasta_path + '.fai', sample_name)

    additional_headers = [
        '##INFO=<ID=CLASS,Number=1,Type=String,Description="Variant class prediction from legacy strict/lenient filters">',
        '##INFO=<ID=CONFIDENCE_SCORE,Number=1,Type=Float,Description="Confidence score based on multiple evidence sources">'
    ]

    header_lines = vcf_header.strip().split('\n')
    insert_pos = len(header_lines) - 1
    for header in additional_headers:
        header_lines.insert(insert_pos, header)
        insert_pos += 1

    return '\n'.join(header_lines) + '\n'


def write_somatic_vcf(passed_breakpoints, output_filename, ref_fasta_path):
    if not passed_breakpoints:
        print("No variants passed filtering.")
        return

    vcf_header = generate_somatic_vcf_header(ref_fasta_path)
    ref_fasta = pysam.FastaFile(ref_fasta_path)

    with open(output_filename, 'w') as f:
        f.write(vcf_header)
        for bp in passed_breakpoints:
            vcf_line = generate_vcf_line(bp, ref_fasta)
            f.write(vcf_line)

    ref_fasta.close()
    print(f"Filtered variants written to: {output_filename}")