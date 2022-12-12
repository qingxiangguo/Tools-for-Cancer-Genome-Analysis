# _*_ coding=utf-8 _*_

with open("/projects/b1171/qgn1237/4_single_cell_SV_chimera/1_smooth_seq_95_sc_K562_SMRT/coverage_list", "a+") as f, \
        open("/home/qgn1237/qgn1237/4_single_cell_SV_chimera/1_smooth_seq_95_sc_K562_SMRT/good_list", "a+") as f2:
    lines_list = f.readlines()
    for line in lines_list:
        line = line.strip()  # You need to remove \n, because the boolean value of \n line is True
        if line:  # You can not delete this if, even if the \n is deleted, you need to get rid of empty list
            if float(line.split()[1]) > 0.10:  # You need to change the type before comparing
                f2.write(line.split()[0] + "\t" + line.split()[1] + "\n")
