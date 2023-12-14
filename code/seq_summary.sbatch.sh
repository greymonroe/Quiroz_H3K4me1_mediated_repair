#!/bin/bash -l
#SBATCH -o /home/gmonroe/slurm-log/%j-stdout.txt
#SBATCH -e /home/gmonroe/slurm-log/%j-stderr.txt
#SBATCH -J seq_sum
#SBATCH -t 96:00:00
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=32G
#SBATCH --partition=bmh
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gmonroe@ucdavis.edu
#!/bin/bash


module load samtools

# Define the output file
output_file="read_counts_depth.tsv"
echo -e "PREFIX\tRaw_Reads_R1\traw_read_count_1" > "$output_file"

# Loop over each PREFIX
for PREFIX in 1B_S1 1C_S2 1D_S3 1E_S4 1F_S5 1G_S6 1H_S7 1I_S8 2C_S9 2D_S10 2E_S11 2F_S12 2I_S13 2L_S14 5B_S15 5I_S16 5J_S17; do
echo $PREFIX
# Count raw reads for R1 and R3
raw_read_count_1=$(zgrep -c '^@' "0_raw_data/Data/ngpqbyc4ib/Un_DTSA603/Project_GMDQ_Nova629P_Quiroz/${PREFIX}_L002_R1_001.fastq.gz")

# Count mapped reads
mapped_read_count=$(samtools view -c -F 4 "2_bam_fix_markdup/$PREFIX.fix.markdup.bam")

# Write the counts and depths to the output file
echo -e "${PREFIX}\t${raw_read_count_1}\t${mapped_read_count}" >> "$output_file"
done

# Print the contents of the file
cat "$output_file"

#!/bin/bash -l
#SBATCH -o /home/gmonroe/slurm-log/%j-stdout.txt
#SBATCH -e /home/gmonroe/slurm-log/%j-stderr.txt
#SBATCH -J seq_sum
#SBATCH -t 96:00:00
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=32G
#SBATCH --partition=bmh
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gmonroe@ucdavis.edu
#!/bin/bash


module load samtools

# Define the output file
output_file="read_counts_depth.tsv"
echo -e "PREFIX\tRaw_Reads_R1\tRaw_Reads_R3\tTrimmed_Reads_R1\tTrimmed_Reads_R2\tMapped_Reads\tUnique_NonDup_Reads\tTotal_Raw_Depth\tTotal_NonDup_Depth" > "$output_file"

# Loop over each PREFIX
for PREFIX in 1B_S1 1C_S2 1D_S3 1E_S4 1F_S5 1G_S6 1H_S7 1I_S8 2C_S9 2D_S10 2E_S11 2F_S12 2I_S13 2L_S14 5B_S15 5I_S16 5J_S17; do
# Count raw reads for R1 and R3
raw_read_count_1=$(zgrep -c '^@' "0_raw_data/Data/ngpqbyc4ib/Un_DTSA603/Project_GMDQ_Nova629P_Quiroz/${PREFIX}_L002_R1_001.fastq.gz")
raw_read_count_3=$(zgrep -c '^@' "0_raw_data/Data/ngpqbyc4ib/Un_DTSA603/Project_GMDQ_Nova629P_Quiroz/${PREFIX}_L002_R3_001.fastq.gz")

# Count trimmed reads for R1 and R2
trimmed_read_count_1=$(zgrep -c '^@' "1_fastq/${PREFIX}_1.trimmed.fastq.gz")
trimmed_read_count_2=$(zgrep -c '^@' "1_fastq/${PREFIX}_2.trimmed.fastq.gz")

# Count mapped reads
mapped_read_count=$(samtools view -c -F 4 "2_bam_fix_markdup/$PREFIX.fix.markdup.bam")

# Count unique mapped non-duplicated reads
unique_non_dup_read_count=$(samtools view -c -F 0x400 -F 0x100 "2_bam_fix_markdup/$PREFIX.fix.markdup.bam")

 # Calculate total raw read depth
total_raw_depth=$(samtools depth -a "2_bam_fix_markdup/$PREFIX.fix.markdup.bam" | awk '{sum+=$3} END {print sum}')

# Calculate total non-duplicated read depth
total_non_dup_depth=$(samtools view -F 0x400 -F 0x100 "2_bam_fix_markdup/$PREFIX.fix.markdup.bam" | \
samtools depth -a - | awk '{sum+=$3} END {print sum}')

# Write the counts and depths to the output file
echo -e "${PREFIX}\t${raw_read_count_1}\t${raw_read_count_3}\t${trimmed_read_count_1}\t${trimmed_read_count_2}\t${mapped_read_count}\t${unique_non_dup_read_count}\t${total_raw_depth}\t${total_non_dup_depth}" >> "$output_file"
done


# Print the contents of the file
cat "$output_file"



#!/bin/bash -l
#SBATCH -o /home/gmonroe/slurm-log/%j-stdout.txt
#SBATCH -e /home/gmonroe/slurm-log/%j-stderr.txt
#SBATCH -J seq_sum
#SBATCH -t 96:00:00
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --partition=bmh
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gmonroe@ucdavis.edu

module load samtools

output_file="read_counts_depth.tsv"
echo -e "PREFIX\tRaw_Reads_R1\tRaw_Reads_R3\tTrimmed_Reads_R1\tTrimmed_Reads_R2\tMapped_Reads\tUnique_NonDup_Reads\tTotal_Raw_Depth\tTotal_NonDup_Depth" > "$output_file"

process_count=0
max_processes=8 # You may adjust this number based on your resource availability

calculate_counts_and_depth() {
  
  PREFIX=$1
  # Count raw reads for R1 and R3
  raw_read_count_1=$(zgrep -c '^@' "0_raw_data/Data/ngpqbyc4ib/Un_DTSA603/Project_GMDQ_Nova629P_Quiroz/${PREFIX}_L002_R1_001.fastq.gz")
  raw_read_count_3=$(zgrep -c '^@' "0_raw_data/Data/ngpqbyc4ib/Un_DTSA603/Project_GMDQ_Nova629P_Quiroz/${PREFIX}_L002_R3_001.fastq.gz")
  
  # Count trimmed reads for R1 and R2
  trimmed_read_count_1=$(zgrep -c '^@' "1_fastq/${PREFIX}_1.trimmed.fastq.gz")
  trimmed_read_count_2=$(zgrep -c '^@' "1_fastq/${PREFIX}_2.trimmed.fastq.gz")
  
  # Count mapped reads
  mapped_read_count=$(samtools view -c -F 4 "2_bam_fix_markdup/$PREFIX.fix.markdup.bam")
  
  # Count unique mapped non-duplicated reads
  unique_non_dup_read_count=$(samtools view -c -F 0x400 -F 0x100 "2_bam_fix_markdup/$PREFIX.fix.markdup.bam")
  
 # Calculate total raw read depth
total_raw_depth=$(samtools depth -a "2_bam_fix_markdup/$PREFIX.fix.markdup.bam" | awk '{sum+=$3} END {print sum}')

# Calculate total non-duplicated read depth
total_non_dup_depth=$(samtools view -F 0x400 -F 0x100 "2_bam_fix_markdup/$PREFIX.fix.markdup.bam" | \
samtools depth -a - | awk '{sum+=$3} END {print sum}')

  echo -e "${PREFIX}\t${raw_read_count_1}\t${raw_read_count_3}\t${trimmed_read_count_1}\t${trimmed_read_count_2}\t${mapped_read_count}\t${unique_non_dup_read_count}\t${total_raw_depth}\t${total_non_dup_depth}" >> "temp_$PREFIX.tsv"
}

export -f calculate_counts_and_depth

for PREFIX in 1B_S1 1C_S2 1D_S3 1E_S4 1F_S5 1G_S6 1H_S7 1I_S8 2C_S9 2D_S10 2E_S11 2F_S12 2I_S13 2L_S14 5B_S15 5I_S16 5J_S17; do
calculate_counts_and_depth $PREFIX &
  let process_count+=1
# If we've reached the maximum number of processes, wait for all to complete
if (( process_count == max_processes )); then
wait
process_count=0
fi
done

# Wait for any remaining background processes to finish
wait

# Now concatenate all the temporary files into the final output file
for PREFIX in 1B_S1 1C_S2 1D_S3 1E_S4 1F_S5 1G_S6 1H_S7 1I_S8 2C_S9 2D_S10 2E_S11 2F_S12 2I_S13 2L_S14 5B_S15 5I_S16 5J_S17; do
cat "temp_$PREFIX.tsv" >> "$output_file"
rm "temp_$PREFIX.tsv"
done

# Print the contents of the file
cat "$output_file"


