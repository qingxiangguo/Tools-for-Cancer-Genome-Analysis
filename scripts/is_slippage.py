#!/usr/bin/python

def is_slippage(chrm, pos, ref, alt, indel_type, genome_seq):
    if indel_type == "INS":
        sequence_to_check = alt[1:]
    elif indel_type == "DEL":
        sequence_to_check = ref[1:]

    if repeat_checker(seq_to_check):
        pattern = rf"({seq_to_check[0]}){{4,}}"
    else:
        pattern = rf"({seq_to_check}){{4,}}"

    pat = re.compile(pattern)
    match = pat.search(genome_seq[chrm][pos + len(sequence_to_check) - 1 : pos + len(sequence_to_check)*4])

    if match:
        return True
    else:
        return False
