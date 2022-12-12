# _*_ coding=utf-8 _*_
import subprocess

# Insert your SRR accession numbers here

sra_accession_number = ["SRR11951487", "SRR11951488", "SRR11951489", "SRR11951490",
                        "SRR11951491", "SRR11951492", "SRR11951493", "SRR11951494", "SRR11951495", "SRR11951496", "SRR11951497",
                        "SRR11951498", "SRR11951499", "SRR11951500", "SRR11951501", "SRR11951502", "SRR11951503", "SRR11951504",
                        "SRR11951505", "SRR11951506", "SRR11951507", "SRR11951508", "SRR11951509", "SRR11951510", "SRR11951511",
                        "SRR11951512", "SRR11951513", "SRR11951514", "SRR11951515", "SRR11951516", "SRR11951517", "SRR11951519",
                        "SRR11951520", "SRR11951521", "SRR11951522", "SRR11951523", "SRR11951524", "SRR11951525", "SRR11951526",
                        "SRR11951527", "SRR11951528", "SRR11951530", "SRR11951531", "SRR11951532", "SRR11951533", "SRR11951534",
                        "SRR11951535", "SRR11951536", "SRR11951537", "SRR11951542"]

for sra_id in sra_accession_number:
        minimap_cmd = "minimap2 -ax map-hifi /home/qgn1237/qgn1237/1_my_database/GRCh38_p13/minimap2_index/GRCh38.p13.genome.mmi /projects/b1171/qgn1237/2_raw_data/smooth_seq_95_sc_K562_SMRT/" + sra_id + "/" + sra_id + ".fastq" + " > /home/qgn1237/qgn1237/4_single_cell_SV_chimera/1_smooth_seq_95_sc_K562_SMRT/" + sra_id + "/" + sra_id + ".sam" + " -t 12"
        subprocess.call(minimap_cmd, shell=True)
