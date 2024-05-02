

module load bwa/0.7.17

# make the common path available in variable (this is optional, and just for the looks
dir=/home/audett/projects/def-idworkin/audett/SSD/genomes/F100/
out=/home/audett/scratch/SSD/F100/mapped

# make it arrays so you can index them
declare -a forward=( /home/audett/projects/def-idworkin/audett/SSD/genomes/F100/*_p1_seq.txt )

# get the individual reads.
R1=${forward[${SLURM_ARRAY_TASK_ID}]}
 
# get the prefix of the reads
#prefix=${R1%%_R1_trimmed.fastq}
 
name=${R1}
base=`basename ${name} _p1_seq.txt`
 
bwa mem -t ${SLURM_CPUS_PER_TASK} -M /home/audett/projects/def-idworkin/audett/SSD/ref/all_ref.fa ${dir}${base}_p1_seq.txt \
   ${dir}${base}_p2_seq.txt > ${out}/${base}.sam

# -M Mark shorter split hits as secondary (for Picard compatibility).
# -t used to tell BWA to use 32 threads to speed things up
