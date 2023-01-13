#!/usr/bin/env fish
for file in (ls *.slurm)
    sbatch $file
end