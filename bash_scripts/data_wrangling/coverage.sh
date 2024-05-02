

module load samtools

parallel --will-cite --jobs 8 samtools depth \
::: ./*.coverage \
> ./{}
