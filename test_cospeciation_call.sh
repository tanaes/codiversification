echo 'Test summarized by host species'

#summarize_otu_by_cat.py -i Example_input/sample_map.txt -c Example_input/otu_table.txt -m HostSpecies -o Example_input/otu_table_HostSpecies.txt
#single_rarefaction.py -i Example_input/otu_table_HostSpecies.txt -o Example_input/otu_table_HostSpecies_rarified.txt -d 100
#filter_otu_table.py -i Example_input/otu_table_HostSpecies_rarified.txt -o Example_input/otu_table_HostSpecies_rarified_filtered.txt -c 3 -s 3

#subcluster OTUs into cOTUs

#python ./otu_subcluster.py \
#-i Example_input/otu_map.txt \
#-o Example_output/new_cOTUs \
#-f Example_input/seqs.fna \
#-p Example_input/99_fasttree_muscle_params.txt \
#-t Example_input/otu_table_rarified_filtered.biom \
#--force

#run cospeciation test on each pOTU

python ./test_cospeciation.py \
-i Example_output/new_cOTUs \
-p Example_input/otu_table_rarified_filtered.biom \
--host_tree_fp Example_input/host_tree_SampleID.tre \
-o Example_output/new_hommola_SampleID \
-T hommola \
-t Example_input/taxonomy.txt \
-m Example_input/sample_map.txt \
--force

#summarize results & run multiple test correction

#python ./summarize_results.py \
#-i Example_output/cOTUs_HostSpecies \
#-r Example_output/hommola_test_HostSpecies \
#-o Example_output/hommola_test_HostSpecies_corrected \
#-t Example_input/taxonomy.txt \
#--force

