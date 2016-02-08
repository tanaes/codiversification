widths=(91 94 97)

# Input base analysis directory info
base_dir=/home/jgsanders/base_dir
collapse_field=tree_names

cd $base_dir


for width in ${widths[@]}
do

cat > ${base_dir}/${width}_analysis.sh << EOF
#!/bin/bash 
#PBS -k oe
#PBS -N ${width}_subcl
#PBS -l nodes=1:ppn=16
#PBS -l mem=16gb
#PBS -l walltime=12:00:00
#PBS -m abe
#PBs -M jonsan@gmail.com

source /home/jgsanders/.bashrc
source activate qiime

width=${width}
collapse_field=${collapse_field}
base_dir=${base_dir}

sample_map=${base_dir}/host_map.txt
host_tree=${base_dir}/hosts.tre
fasta_input=${base_dir}/seqs.fna
scripts_dir=/home/jgsanders/git_sw/codiversification
scripts_bin=/home/jgsanders/git_sw/script_bin
subcluster_params=/home/jgsanders/git_sw/codiversification/params/99_fasttree_muscle_sumaclust_params.txt


EOF

cat >> ${base_dir}/${width}_analysis.sh  << 'EOF'

output_dir=${base_dir}/${width}

# De-novo cluster at different thresholds
# methods: cd-hit, swarm, usearch61

mkdir ${output_dir}
cd ${output_dir}

#### run collapsing by HostSpecies ####

# collapse by Species
collapse_samples.py -b ${base_dir}/host_${width}_otu_table.biom \
-m ${sample_map} \
--output_biom_fp ${output_dir}/otu_table.${collapse_field}.biom \
--output_mapping_fp ${output_dir}/map.${collapse_field}.txt \
--collapse_fields ${collapse_field}

# filter otu table for the host species we want
# eval filter_samples_from_otu_table.py -i ${output_dir}/otu_table_${collapse_field}.biom \
# -o ${output_dir}/otu_table_${collapse_field}_${filter_field}-${filter_value}.biom \
# -m ${output_dir}/map_${collapse_field}.txt \
# -s '${filter_field}:${filter_value}'

# rarify for even sampling across species
single_rarefaction.py -i ${output_dir}/otu_table.${collapse_field}.biom \
-o ${output_dir}/otu_table.${collapse_field}.10000.biom \
-d 10000

# subset filtered otu table by abundance
filter_otus_from_otu_table.py \
-i ${output_dir}/otu_table.${collapse_field}.10000.biom \
-s 3 \
-o ${output_dir}/otu_table.${collapse_field}.10000.s3.biom

# split OTU table for parallel processing
mkdir ${output_dir}/temp_biom
cd ${output_dir}/temp_biom
python ${scripts_bin}/split_biom.py -i ${output_dir}/otu_table.${collapse_field}.10000.s3.biom -n 16
cd ${output_dir}

mkdir ${output_dir}/subclustered_otus

# OTU subcluster
eval parallel 'python ${scripts_dir}/otu_subcluster.py \
-i ${base_dir}/final_otu_map_mc2_${width}.txt \
-o ${output_dir}/subclustered_otus \
-f ${fasta_input} \
-p ${subcluster_params} \
--force \
-b {}' ::: ${output_dir}/temp_biom/chunk*.biom


# python ${scripts_dir}/otu_subcluster.py \
# -i ${base_dir}/final_otu_map_mc2_${width}.txt \
# -o ${output_dir}/subclustered_otus \
# -f ${fasta_input} \
# -p ${subcluster_params} \
# --force \
# -b ${output_dir}/otu_table.${collapse_field}.10000.s3.biom

#run cospeciation test
echo "source ~/.bashrc 
source activate qiime
python ${scripts_dir}/test_cospeciation.py \
-i ${output_dir}/subclustered_otus \
-p ${output_dir}/otu_table.${collapse_field}.10000.s3.biom \
--host_tree_fp ${host_tree} \
-o ${output_dir}/hommola_test_${collapse_field} \
-T hommola_host \
-m ${sample_map} \
--collapse_fields ${collapse_field} \
--permutations 10000 \
--force" | \
qsub -k oe \
-N ${width}_cosp \
-l nodes=1:ppn=1 \
-l pmem=8gb \
-l walltime=24:00:00 \
-m abe \
-M jonsan@gmail.com

EOF

chmod 755 ${base_dir}/${width}_analysis.sh

done





