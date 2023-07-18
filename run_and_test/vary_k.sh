
for sub_dir in spring minicom pgrc; do {
	cd /vol1/agis/ruanjue_group/jifahu/zsmem-experiments/${sub_dir} || exit
	pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/bwamem_seeding -k15 -t16 ../bwamem/index_celegans/ref elegans_SRR16905161.reads 2>&1 | tee seed_log/elegans.k15.bwa
	pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/bwamem_seeding -k17 -t16 ../bwamem/index_celegans/ref elegans_SRR16905161.reads 2>&1 | tee seed_log/elegans.k17.bwa
	pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/bwamem_seeding -k19 -t16 ../bwamem/index_celegans/ref elegans_SRR16905161.reads 2>&1 | tee seed_log/elegans.k19.bwa
	pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/bwamem_seeding -k21 -t16 ../bwamem/index_celegans/ref elegans_SRR16905161.reads 2>&1 | tee seed_log/elegans.k21.bwa
	pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/bwamem_seeding -k23 -t16 ../bwamem/index_celegans/ref elegans_SRR16905161.reads 2>&1 | tee seed_log/elegans.k23.bwa
} done
wait
pstat -m "Vary-K-BWA"

for sub_dir in spring minicom pgrc; do {
	cd /vol1/agis/ruanjue_group/jifahu/zsmem-experiments/${sub_dir} || exit
	pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/comp_seeding -k15 -t16 ../bwamem/index_celegans/ref elegans_SRR16905161.reads 2>&1 | tee seed_log/elegans.k15.comp
	pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/comp_seeding -k17 -t16 ../bwamem/index_celegans/ref elegans_SRR16905161.reads 2>&1 | tee seed_log/elegans.k17.comp
	pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/comp_seeding -k19 -t16 ../bwamem/index_celegans/ref elegans_SRR16905161.reads 2>&1 | tee seed_log/elegans.k19.comp
	pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/comp_seeding -k21 -t16 ../bwamem/index_celegans/ref elegans_SRR16905161.reads 2>&1 | tee seed_log/elegans.k21.comp
	pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/comp_seeding -k23 -t16 ../bwamem/index_celegans/ref elegans_SRR16905161.reads 2>&1 | tee seed_log/elegans.k23.comp
} done
wait
pstat -m "Vary-K-COMP"

