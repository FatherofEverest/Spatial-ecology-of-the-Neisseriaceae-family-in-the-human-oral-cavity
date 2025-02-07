pwd=/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter
python SCRIPTS/p-get_percent_identity.py \
       -B $pwd/28_MAPPING_QUALITY/soi.txt \
       -p \
       -o $pwd/28_MAPPING_QUALITY/PERCENT_IDENTITY_CORE_GENES_H_parainfluenzae_str_M1C160_1_id_GCA_014931275_1.txt \
       -range 65,100 \
       -autobin \
       -m \
       -melted \
       -G goi \
       -a $pwd/11_GENE_CALLS/gene_calls_summary.txt \
       -interpolate 70,100,400 \
       -x

