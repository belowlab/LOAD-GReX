gwas_parsing.py \
        -gwas_file $1'_filtered.txt' \
        -liftover ~/tools/MASHR/data/liftover/hg19ToHg38.over.chain.gz \
        -snp_reference_metadata ~/tools/MASHR/data/reference_panel_1000G/variant_metadata.txt.gz METADATA \
        -output_column_map rsID variant_id \
        -output_column_map Reference_allele non_effect_allele \
        -output_column_map Effect_allele effect_allele \
        -output_column_map BETA effect_size \
        -output_column_map P pvalue \
        -output_column_map Chr chromosome \
        --chromosome_format \
        -output_column_map Position position \
        -output_column_map EAF_U frequency \
        -output_column_map SE standard_error \
        --insert_value sample_size $2 --insert_value n_cases $3 \
        -output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases\
        -output $1.harmonized.txt.gz

SPrediXcan.py \
  --model_db_path ~/tools/MASHR/eqtl/mashr/$db.db \
  --gwas_file $in.harmonized.txt.gz \
  --snp_column panel_variant_id \
  --effect_allele_column effect_allele \
  --non_effect_allele_column non_effect_allele \
  --pvalue_column pvalue \ 
  --beta_column effect_size \ 
  --covariance ~/tools/MASHR/eqtl/mashr/$db.txt.gz \ 
  --keep_non_rsid --additional_output --model_db_snp_key varID --throw \
  --output_file $out
