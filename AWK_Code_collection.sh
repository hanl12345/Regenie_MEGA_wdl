#Convert Ted's phecode csv file to regenie input file with Tab as splitter,convert missing value as NA
awk 'NR==FNR {keys[$1]; next} FNR==1 {OFS="\t"; $1=$1; print $1,$0; next} $1 in keys {for (i=1; i<=NF; i++) if ($i == "") $i="NA"; OFS="\t"; $1=$1; print $1,$0}' ./AGD35K_MAC100MR005QCed_Autosome.fam FS="," ./AGD_35k_code_data_floor_wide_format_phecode_x.csv >./AGD_35k_code_data_floor_wide_format_phecode_x_gwas.txt

