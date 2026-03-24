#!/bin/bash
#SBATCH -J compare_v3
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -o compare_v3_%j.out

python compare_got_and_want_v3.py experiments.bk/exclaim_ape_R02B04_dt8_g0008_R02B05 -o dt8_g0008_R02B05.db -j 32
