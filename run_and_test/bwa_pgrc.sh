cd /vol1/agis/ruanjue_group/jifahu/zsmem-experiments/pgrc || exit
pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/bwamem_seeding -t16 ../bwamem/index_ecoli/ref ecoli_SRR1562082_1.reads 2>&1 | tee seed_log/ecoli_SRR1562082_1.bwa
pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/bwamem_seeding -t16 ../bwamem/index_celegans/ref elegans_SRR16905161.reads 2>&1 | tee seed_log/elegans_SRR16905161.bwa
pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/bwamem_seeding -t16 ../bwamem/index_chicken/ref chicken_SRR13537343_1.reads 2>&1 | tee seed_log/chicken_SRR13537343_1.bwa
pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/bwamem_seeding -t16 ../bwamem/index_HG19/ref human_ERP001775_s1_1.reads 2>&1 | tee seed_log/human_ERP001775_s1_1.bwa
pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/bwamem_seeding -t16 ../bwamem/index_HG19/ref human_ERP001775_s2_1.reads 2>&1 | tee seed_log/human_ERP001775_s2_1.bwa
pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/bwamem_seeding -t16 ../bwamem/index_HG19/ref human_ERP001775_s3_1.reads 2>&1 | tee seed_log/human_ERP001775_s3_1.bwa
pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/bwamem_seeding -t16 ../bwamem/index_HG19/ref human_ERP001775_s4_1.reads 2>&1 | tee seed_log/human_ERP001775_s4_1.bwa
pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/bwamem_seeding -t16 ../bwamem/index_HG19/ref human_ERP001775_s5_1.reads 2>&1 | tee seed_log/human_ERP001775_s5_1.bwa
pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/bwamem_seeding -t16 ../bwamem/index_HG19/ref human_ERR194146_1.reads 2>&1 | tee seed_log/human_ERR194146_1.bwa
pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/bwamem_seeding -t16 ../bwamem/index_HG19/ref human_ERR194161_1.reads 2>&1 | tee seed_log/human_ERR194161_1.bwa
pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/bwamem_seeding -t16 ../bwamem/index_HG19/ref human_ERR3239279_1.reads 2>&1 | tee seed_log/human_ERR3239279_1.bwa
pstat /vol1/agis/ruanjue_group/jifahu/project/comp_seed/build/bwamem_seeding -t16 ../bwamem/index_HG19/ref human_SRR10965089_1.reads 2>&1 | tee seed_log/human_SRR10965089_1.bwa
