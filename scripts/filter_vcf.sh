#!/bin/bash

threshold=$1
input_vcf=$2
prefix="${input_vcf%.var.vcf}" # Delete the var.vcf in input
output_vcf="${prefix}.filtered.var.vcf"

bcftools filter -i "FILTER == \"PASS\" && FORMAT/AD[0:1] >= ${threshold}" "${input_vcf}" > "${output_vcf}"