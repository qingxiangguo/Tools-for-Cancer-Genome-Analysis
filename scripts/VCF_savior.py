#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys
import re
import logging
from datetime import datetime
import gzip

class VCFSavior:
    def __init__(self, input_vcf, output_vcf, genome_version=None):
        """
        Initialize VCF Savior
        
        Args:
            input_vcf (str): Input VCF file path
            output_vcf (str): Output VCF file path
            genome_version (str): Genome version (37 or 38)
        """
        self.input_vcf = input_vcf
        self.output_vcf = output_vcf
        self.genome_version = genome_version
        self.temp_files = []
        
        # Setup logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        self.logger = logging.getLogger(__name__)
        
    def cleanup(self):
        """Remove temporary files"""
        for temp_file in self.temp_files:
            if os.path.exists(temp_file):
                os.remove(temp_file)
                
    def create_temp_filename(self, suffix):
        """Create temporary filename"""
        temp_file = f"{self.input_vcf}.temp_{suffix}"
        self.temp_files.append(temp_file)
        return temp_file

    def standardize_chrom_name(self, chrom):
        """
        Standardize chromosome names based on genome version
        Only standardize regular chromosomes (1-22, X, Y, M)
        Leave other contigs unchanged
        
        Args:
            chrom (str): Chromosome name
        Returns:
            str: Standardized chromosome name
        """
        if not self.genome_version:
            return chrom
            
        # Define standard chromosomes
        std_chromosomes = set([str(i) for i in range(1, 23)] + ['X', 'Y', 'M'])
        
        # For GRCh38, add 'chr' prefix only to standard chromosomes
        if self.genome_version == '38':
            if chrom.startswith('chr'):
                stripped_chrom = chrom[3:]
                # Only keep 'chr' prefix if it's a standard chromosome
                if stripped_chrom in std_chromosomes:
                    return chrom
                return stripped_chrom
            else:
                # Add 'chr' only if it's a standard chromosome
                if chrom in std_chromosomes:
                    return f"chr{chrom}"
                return chrom
                
        # For GRCh37, remove 'chr' prefix only from standard chromosomes
        elif self.genome_version == '37':
            if chrom.startswith('chr'):
                stripped_chrom = chrom[3:]
                # Only remove 'chr' prefix if it's a standard chromosome
                if stripped_chrom in std_chromosomes:
                    return stripped_chrom
                return chrom
            return chrom
            
        return chrom

    def process_header_line(self, line):
        """
        Process header line for chromosome standardization
        
        Args:
            line (str): Header line
        Returns:
            str: Processed header line
        """
        if not self.genome_version:
            return line
            
        if line.startswith('##contig=<ID='):
            chrom = re.search(r'ID=([^,>]+)', line).group(1)
            std_chrom = self.standardize_chrom_name(chrom)
            return line.replace(f'ID={chrom}', f'ID={std_chrom}')
        return line

    def process_vcf_line(self, line):
        """
        Process VCF data line for chromosome standardization
        
        Args:
            line (str): VCF data line
        Returns:
            str: Processed VCF data line
        """
        if not self.genome_version or line.startswith('#'):
            return line
            
        fields = line.strip().split('\t')
        if len(fields) > 0:
            # Standardize main chromosome
            fields[0] = self.standardize_chrom_name(fields[0])
            
            # Check INFO field for CHR2 or other chromosome references
            if len(fields) > 7:
                info_fields = fields[7].split(';')
                for i, info in enumerate(info_fields):
                    if info.startswith(('CHR2=', 'CHROM2=')):
                        key, value = info.split('=')
                        std_value = self.standardize_chrom_name(value)
                        info_fields[i] = f'{key}={std_value}'
                fields[7] = ';'.join(info_fields)
                
        return '\t'.join(fields)

    def extract_header_definitions(self, header_lines):
        """
        Extract existing INFO and FORMAT field definitions
        
        Args:
            header_lines (list): List of header lines
        Returns:
            tuple: Sets of defined INFO and FORMAT fields
        """
        info_fields = set()
        format_fields = set()
        
        for line in header_lines:
            if line.startswith('##INFO=<ID='):
                match = re.search(r'ID=([^,]+)', line)
                if match:
                    info_fields.add(match.group(1))
            elif line.startswith('##FORMAT=<ID='):
                match = re.search(r'ID=([^,]+)', line)
                if match:
                    format_fields.add(match.group(1))
        
        return info_fields, format_fields

    def find_undefined_fields(self, vcf_lines, defined_info, defined_format):
        """
        Find undefined INFO and FORMAT fields
        
        Args:
            vcf_lines (list): List of VCF data lines
            defined_info (set): Set of defined INFO fields
            defined_format (set): Set of defined FORMAT fields
        Returns:
            tuple: Sets of undefined INFO and FORMAT fields
        """
        undefined_info = set()
        undefined_format = set()
        
        for line in vcf_lines:
            if isinstance(line, str):
                fields = line.strip().split('\t')
            else:
                fields = line
                
            # Check INFO fields
            info = fields[7]
            if info != '.':
                for item in info.split(';'):
                    if '=' in item:
                        field = item.split('=')[0]
                        if field not in defined_info:
                            undefined_info.add(field)
                    elif item not in defined_info:
                        undefined_info.add(item)
            
            # Check FORMAT fields
            if len(fields) > 8:
                format_fields = fields[8].split(':')
                for field in format_fields:
                    if field not in defined_format:
                        undefined_format.add(field)
        
        return undefined_info, undefined_format

    def fix_header_definitions(self, header_lines):
        """
        Fix header definitions
        
        Args:
            header_lines (list): List of header lines
        Returns:
            list: Fixed header lines
        """
        fixed_headers = []
        pl_defined = False
        
        for line in header_lines:
            if line.startswith('##FORMAT=<ID=PL,'):
                line = '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes">'
                pl_defined = True
            elif line.startswith('##FORMAT=<ID=GT,') and any(h.startswith('##FORMAT=<ID=GT,') for h in fixed_headers):
                continue
                
            fixed_headers.append(line)
        
        if not pl_defined and any(h.startswith('##FORMAT=<ID=') for h in fixed_headers):
            fixed_headers.insert(-1, '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes">')
        
        return fixed_headers

    def fix_genotype_field(self, format_str, sample_str):
        """
        Fix genotype field ensuring valid format
        
        Args:
            format_str (str): FORMAT column string
            sample_str (str): Sample column string
        Returns:
            tuple: Fixed FORMAT and sample strings
        """
        format_fields = format_str.split(':')
        sample_fields = sample_str.split(':')
        
        if 'GT' not in format_fields:
            format_fields.insert(0, 'GT')
            sample_fields.insert(0, '1/1')
        elif format_fields[0] != 'GT':
            gt_idx = format_fields.index('GT')
            format_fields.insert(0, format_fields.pop(gt_idx))
            sample_fields.insert(0, sample_fields.pop(gt_idx))
        
        if sample_fields[0] in ['./.', '.', './1', '1/.']:
            sample_fields[0] = '1/1'
        
        return ':'.join(format_fields), ':'.join(sample_fields)

    def fix_svlen(self, info_field, pos):
        """
        Fix missing or invalid SVLEN values in INFO field
        
        Args:
            info_field (str): INFO field from VCF
            pos (str): Position field from VCF
        Returns:
            str: Fixed INFO field
        """
        info_dict = {}
        for kv in info_field.split(';'):
            if '=' in kv:
                k, v = kv.split('=', 1)
                info_dict[k] = v
            else:
                info_dict[kv] = True

        svlen = info_dict.get("SVLEN", None)
        svtype = info_dict.get("SVTYPE", None)
        curr_pos = int(pos)
        end_str = info_dict.get("END", None)

        if svlen is None or svlen == ".":
            new_svlen = None
            if end_str is not None and end_str.isdigit():
                end = int(end_str)
                if svtype == "DEL":
                    new_svlen = -(end - curr_pos + 1)
                elif svtype == "INS":
                    calc_len = end - curr_pos
                    if calc_len < 1:
                        calc_len = 1
                    new_svlen = calc_len
                elif svtype in ["DUP", "INV"]:
                    calc_len = end - curr_pos + 1
                    if calc_len < 1:
                        calc_len = 1
                    new_svlen = calc_len
                else:
                    new_svlen = 0
            else:
                if svtype == "INS":
                    new_svlen = 1
                elif svtype == "DEL":
                    new_svlen = -1
                else:
                    new_svlen = 0
            
            info_dict["SVLEN"] = str(new_svlen)

        new_info = []
        for k, v in info_dict.items():
            if v is True:
                new_info.append(k)
            else:
                new_info.append(f"{k}={v}")
        
        return ';'.join(new_info)

    def fix_vcf(self):
        """Main process to fix VCF file"""
        try:
            self.logger.info(f"Starting VCF processing with genome version {self.genome_version}")
            
            opener = gzip.open if self.input_vcf.endswith('.gz') else open
            mode = 'rt' if self.input_vcf.endswith('.gz') else 'r'
            
            with opener(self.input_vcf, mode) as f:
                content = f.readlines()

            header_lines = []
            data_lines = []
            changes_made = set()

            for line in content:
                if line.startswith('#'):
                    processed_line = self.process_header_line(line.strip())
                    if processed_line != line.strip():
                        changes_made.add('chromosome_names_in_header')
                    header_lines.append(processed_line)
                else:
                    processed_line = self.process_vcf_line(line)
                    if processed_line != line:
                        changes_made.add('chromosome_names_in_variants')
                    data_lines.append(processed_line)

            defined_info, defined_format = self.extract_header_definitions(header_lines)
            undefined_info, undefined_format = self.find_undefined_fields(data_lines, defined_info, defined_format)

            if undefined_info:
                changes_made.add('added_missing_info_definitions')
            if undefined_format:
                changes_made.add('added_missing_format_definitions')

            new_headers = []
            for field in sorted(undefined_info):
                new_headers.append(f'##INFO=<ID={field},Number=.,Type=String,Description="Auto-generated definition for {field}">')
            for field in sorted(undefined_format):
                new_headers.append(f'##FORMAT=<ID={field},Number=.,Type=String,Description="Auto-generated definition for {field}">')

            fixed_headers = self.fix_header_definitions(header_lines[:-1] + new_headers)

            with open(self.output_vcf, 'w') as f:
                for header in fixed_headers:
                    f.write(header + '\n')
                f.write(header_lines[-1] + '\n')

                for line in data_lines:
                    fields = line.strip().split('\t')
                    if fields[6] != 'PASS':
                        changes_made.add('filters_set_to_pass')
                    fields[6] = 'PASS'
                    
                    old_info = fields[7]
                    fields[7] = self.fix_svlen(fields[7], fields[1])
                    if old_info != fields[7]:
                        changes_made.add('fixed_svlen_values')
                    
                    if len(fields) >= 10:
                        old_format = fields[8]
                        old_sample = fields[9]
                        fields[8], fields[9] = self.fix_genotype_field(fields[8], fields[9])
                        if old_format != fields[8] or old_sample != fields[9]:
                            changes_made.add('fixed_genotype_fields')
                    
                    f.write('\t'.join(fields) + '\n')

            if changes_made:
                self.logger.info("Changes made to the VCF file:")
                for change in sorted(changes_made):
                    self.logger.info(f"- {change.replace('_', ' ').capitalize()}")
            else:
                self.logger.info("No changes were necessary for this VCF file")

            sorted_vcf = self.output_vcf.replace('.vcf', '_sorted.vcf')
            
            self.logger.info("Sorting VCF file...")
            sort_cmd = f"cat {self.output_vcf} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1V -k2,2n\"}}' > {sorted_vcf}"
            subprocess.run(sort_cmd, shell=True, check=True)

            self.logger.info("Compressing sorted VCF...")
            subprocess.run(f"bgzip -c {sorted_vcf} > {sorted_vcf}.gz", shell=True, check=True)

            self.logger.info("Creating index...")
            subprocess.run(f"tabix -p vcf {sorted_vcf}.gz", shell=True, check=True)

            self.logger.info(f"Processing complete. Files generated:")
            self.logger.info(f"1. Fixed VCF: {self.output_vcf}")
            self.logger.info(f"2. Sorted VCF: {sorted_vcf}")
            self.logger.info(f"3. Compressed VCF: {sorted_vcf}.gz")
            self.logger.info(f"4. Index file: {sorted_vcf}.gz.tbi")

        except Exception as e:
            self.logger.error(f"Error: {str(e)}")
            raise
        finally:
            self.cleanup()

def main():
    parser = argparse.ArgumentParser(description='VCF Savior - Comprehensive VCF fixing tool')
    parser.add_argument('-i', '--input', required=True, help='Input VCF file')
    parser.add_argument('-o', '--output', required=True, help='Output VCF file')
    parser.add_argument('-g', '--genome', choices=['37', '38'], help='Genome version (37 or 38) for chromosome naming')
    
    args = parser.parse_args()
    
    savior = VCFSavior(args.input, args.output, args.genome)
    savior.fix_vcf()

if __name__ == '__main__':
    main()

