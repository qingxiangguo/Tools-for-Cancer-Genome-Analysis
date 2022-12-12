# _*_ coding=utf-8 _*_

import subprocess

# Insert your directories here
sra_accession_number = ["SRR11951439", "SRR11951443", "SRR11951444", "SRR11951445", "SRR11951446", "SRR11951447", "SRR11951448",
                        "SRR11951449", "SRR11951450", "SRR11951451", "SRR11951452", "SRR11951453", "SRR11951454", "SRR11951455",
                        "SRR11951456", "SRR11951457", "SRR11951458", "SRR11951459", "SRR11951460", "SRR11951461", "SRR11951462",
                        "SRR11951463", "SRR11951464", "SRR11951465", "SRR11951466", "SRR11951467", "SRR11951468", "SRR11951469",
                        "SRR11951470", "SRR11951471", "SRR11951472", "SRR11951473", "SRR11951474", "SRR11951475", "SRR11951476",
                        "SRR11951477", "SRR11951478", "SRR11951479", "SRR11951480", "SRR11951481", "SRR11951482", "SRR11951483",
                        "SRR11951484", "SRR11951485", "SRR11951486", "SRR11951487", "SRR11951488", "SRR11951489", "SRR11951490",
                        "SRR11951491", "SRR11951492", "SRR11951493", "SRR11951494", "SRR11951495", "SRR11951496", "SRR11951497",
                        "SRR11951498", "SRR11951499", "SRR11951500", "SRR11951501", "SRR11951502", "SRR11951503", "SRR11951504",
                        "SRR11951505", "SRR11951506", "SRR11951507", "SRR11951508", "SRR11951509", "SRR11951510", "SRR11951511",
                        "SRR11951512", "SRR11951513", "SRR11951514", "SRR11951515", "SRR11951516", "SRR11951517", "SRR11951519",
                        "SRR11951520", "SRR11951521", "SRR11951522", "SRR11951523", "SRR11951524", "SRR11951525", "SRR11951526",
                        "SRR11951527", "SRR11951528", "SRR11951530", "SRR11951531", "SRR11951532", "SRR11951533", "SRR11951534",
                        "SRR11951535", "SRR11951536", "SRR11951537", "SRR11951542"]

for sra_id in sra_accession_number:
    samtools_cmd = "samtools coverage /home/qgn1237/qgn1237/4_single_cell_SV_chimera/1_smooth_seq_95_sc_K562_SMRT/" + sra_id + "/" + sra_id + "_sort.bam | awk 'NR>1 {suma+=$5; sumb+=$3 } END { print suma/sumb}'"  # Sum the 3th column and the 5th column respectively
    output = subprocess.check_output([samtools_cmd], shell=True) # The output of subproces stout is byte style, you need to convert it to str
    output = output.decode("utf-8")
    with open("/home/qgn1237/qgn1237/4_single_cell_SV_chimera/1_smooth_seq_95_sc_K562_SMRT/coverage_list", "a+") as f:
        f.write(sra_id + "\t" + output + "\n")
