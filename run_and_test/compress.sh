
sub_dir=$1
data_id=$2

if [ ! -f /vol1/agis/ruanjue_group/jifahu/dataset/${sub_dir}/${data_id}.fq ]; then
	echo -e 'file does not exist'
	exit
fi

cd /vol1/agis/ruanjue_group/jifahu/zsmem-experiments/spring || exist

pstat /vol1/agis/ruanjue_group/jifahu/software/SPRING/build/spring \
	-c -r -t 16 --no-ids --no-quality \
	-i /vol1/agis/ruanjue_group/jifahu/dataset/${sub_dir}/${data_id}.fq \
	-o archive/${data_id}.spring \
	2>&1 | tee com_log/${data_id}.log

pstat /vol1/agis/ruanjue_group/jifahu/software/SPRING/build/spring \
	-d -t 16 \
	-i archive/${data_id}.spring \
	-o ${data_id}.fq \
	2>&1 | tee dec_log/${data_id}.log

cat ${data_id}.fq | awk 'NR%2==0' > ${sub_dir}_${data_id}.reads
rm ${data_id}.fq