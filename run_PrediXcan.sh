PrediXcan.py \
  --model_db_path ../predixcan_models/gtex_v8/mashr/eqtl/mashr/(*MASHR).db \
  --text_genotypes ../ADGC/chr*.dosage.gz \
  --text_sample_ids ../ADGC/(*study).fam \
  --prediction_output (*MASHR).txt \
  --prediction_summary_output (*MASHR).sum.txt
