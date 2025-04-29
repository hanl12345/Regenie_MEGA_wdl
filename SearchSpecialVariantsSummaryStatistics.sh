tempfile="search_result.txt"
> "$tempfile"
for file in $(gsutil ls gs://fc-secure-540f27be-97ea-4ffd-adb7-c195458eb278/AGD35K_GWAS_SummaryStatistics_0422_2025/AFRRegenie/AGD35kAFR_*.regenie.gz); do
 filename=$(basename "$file")
 # Check if the file contains the target variant (chr1:782105)
 if gsutil cat "$file" | gunzip -c | grep -q "chr1:782105"; then
   echo "Found in $filename:" >> "$tempfile"
   gsutil cat "$file" | gunzip -c | grep "chr1:782105" >> "$tempfile"
 fi
done
