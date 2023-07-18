
for sub_dir in spring minicom pgrc; do {
	cd /vol1/agis/ruanjue_group/jifahu/zsmem-experiments/${sub_dir} || exit
	pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/bwamem_seeding -r0.5 -t16 ../bwamem/index_celegans/ref elegans_SRR16905161.reads 2>&1 | tee seed_log/elegans.r05.bwa
	pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/bwamem_seeding -r1.0 -t16 ../bwamem/index_celegans/ref elegans_SRR16905161.reads 2>&1 | tee seed_log/elegans.r10.bwa
	pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/bwamem_seeding -r1.5 -t16 ../bwamem/index_celegans/ref elegans_SRR16905161.reads 2>&1 | tee seed_log/elegans.r15.bwa
	pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/bwamem_seeding -r2.0 -t16 ../bwamem/index_celegans/ref elegans_SRR16905161.reads 2>&1 | tee seed_log/elegans.r20.bwa
	pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/bwamem_seeding -r2.5 -t16 ../bwamem/index_celegans/ref elegans_SRR16905161.reads 2>&1 | tee seed_log/elegans.r25.bwa
} done
wait
pstat -m "Vary-R-BWA"

for sub_dir in spring minicom pgrc; do {
	cd /vol1/agis/ruanjue_group/jifahu/zsmem-experiments/${sub_dir} || exit
	pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/comp_seeding -r0.5 -t16 ../bwamem/index_celegans/ref elegans_SRR16905161.reads 2>&1 | tee seed_log/elegans.r05.comp
	pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/comp_seeding -r1.0 -t16 ../bwamem/index_celegans/ref elegans_SRR16905161.reads 2>&1 | tee seed_log/elegans.r10.comp
	pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/comp_seeding -r1.5 -t16 ../bwamem/index_celegans/ref elegans_SRR16905161.reads 2>&1 | tee seed_log/elegans.r15.comp
	pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/comp_seeding -r2.0 -t16 ../bwamem/index_celegans/ref elegans_SRR16905161.reads 2>&1 | tee seed_log/elegans.r20.comp
	pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/comp_seeding -r2.5 -t16 ../bwamem/index_celegans/ref elegans_SRR16905161.reads 2>&1 | tee seed_log/elegans.r25.comp
} done
wait
pstat -m "Vary-R-COMP"

