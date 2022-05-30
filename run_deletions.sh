#! /bin/bash
function usage() {
    cat <<USAGE

    Usage: $0 [--contig name] [--reference reference.fa] [--samtools path-to-samtools]
              [--bwa path-to-bwa] [--bedtools path-to-bedtools] [--scripts path-to-scripts-folder]
              [--filepaths path-to-bamfiles]

    Options:
        --contig     name of the contig to be analysed (e.g. Chr1)
        --reference  path to fasta reference used for mapping (e.g. /home/reference.fa)
        --samtools   full path to samtools command (e.g. /home/samtools/samtools)
        --bwa        full path to bwa command (e.g. /home/bwa/bwa)
        --bedtools   full path to bedtools command (e.g. /home/bedtools2/bedtools)
        --scripts    path to folder containing the python scripts (e.g. /home/aDNA-deletions/scripts)
        --filepaths  path to txt file containing sample information
USAGE
    exit 1
}

if [ $# -ne 14 ]; then
    usage
    exit 1
fi

while [ "$1" != "" ]; do
    case $1 in
    --contig)
        shift
        contig_name=$1
        ;;
    --reference)
        shift
        reference_file=$1
        ;;
    --bwa)
        shift
        bwa=$1
        ;;
    --samtools)
        shift
        samtools=$1
        ;;
    --bedtools)
        shift
        bedtools=$1
        ;;
    --scripts)
        shift
        scripts=$1
        ;;
    --filepaths)
        shift
        filepaths=$1
        ;;
    -h | --help)
        usage
        ;;
    *)
        usage
        exit 1
        ;;
    esac
    shift
done

echo "Starting run for contig: "$contig_name

threads=1 #number of threads to use (multi-threading is only used in the first steps of the pipeline)
readlength=35 #minimal readlength used for the ancient genomes (minimal of 35bp is recommended)

$samtools faidx $reference_file
$samtools faidx $reference_file $contig_name > $contig_name".fa"

python $scripts"fasta_to_reads.py" --fasta $contig_name".fa" --readlength $readlength --out $contig_name --target $contig_name

$bwa index $contig_name".fa"
$bwa aln -l 16500 -n 0.01 -o 2 -t $threads $contig_name".fa" $contig_name".fastq.gz" > $contig_name".sai"
$bwa samse $contig_name".fa" $contig_name".sai" $contig_name".fastq.gz" | samtools view -@ $threads -h -Sb - > $contig_name".bam"
rm $contig_name".sai"
$samtools sort -@ $threads $contig_name".bam" > "sorted."$contig_name".bam"
$samtools index -@ $threads "sorted."$contig_name".bam"
rm $contig_name".bam"

echo reference$'\t'sorted.$contig_name.bam > $contig_name".samples.txt"
cat $filepaths >> $contig_name".samples.txt"
sed 1d $contig_name".samples.txt" | cut -f 2 > $contig_name".depth_files.txt"

$samtools depth -aa -Q 30 -r $contig_name "sorted."$contig_name".bam" > "reference."$contig_name".depth"
$samtools depth -aa -Q 30 -r $contig_name -f $contig_name".depth_files.txt" | cut -f 3- > "outgroup."$contig_name".depth"
paste "reference."$contig_name".depth" "outgroup."$contig_name".depth" > $contig_name".depth"

rm "outgroup."$contig_name".depth"
rm "reference."$contig_name".depth"
rm "sorted."$contig_name".bam"
rm "sorted."$contig_name".bam.bai"
rm $contig_name".fa"
rm $contig_name".fa.amb"
rm $contig_name".fa.ann"
rm $contig_name".fa.bwt"
rm $contig_name".fa.pac"
rm $contig_name".fa.sa"

python $scripts"/detect_deletions.py" --depthfile $contig_name".depth" --samples $contig_name".samples.txt"

$bedtools merge -i $contig_name".reference.bed" -d 250 > $contig_name".merged.reference.bed"
$bedtools merge -i $contig_name".outgroup.bed" -d 250 > $contig_name".merged.outgroup.bed"
$bedtools merge -i $contig_name".target.bed" -d 250 > $contig_name".merged.target.bed"

rm $contig_name".reference.bed"
rm $contig_name".outgroup.bed"
rm $contig_name".target.bed"

python $scripts/"merge_bed_files.py" --target $contig_name

$bedtools merge -i $contig_name".intermediate.reference_outgroup.bed" -d 250 | python ~/Desktop/Mammoth_deletions/filter_bed_by_length.py - > $contig_name".reference_outgroup_filter.bed"
$bedtools merge -i $contig_name".intermediate.target.bed" -d 250 | python ~/Desktop/Mammoth_deletions/filter_bed_by_length.py - > $contig_name".target_deletions.bed"

rm $contig_name".intermediate.reference_outgroup.bed"
rm $contig_name".intermediate.target.bed"
rm $contig_name".merged.outgroup.bed"
rm $contig_name".merged.reference.bed"
rm $contig_name".merged.target.bed"
rm $contig_name".depth_files.txt"
rm $contig_name".samples.txt"
rm $contig_name".fastq.gz"

mkdir PDFS

python $scripts"/plot_deletions.py" --depth_file "windowed."$contig_name".depth" \
  --average_coverage $contig_name.average-coverage \
  --outgroup_bed $contig_name".reference_outgroup_filter.bed" \
  --target_bed $contig_name".target_deletions.bed"

mv *$contig_name".window-"*".pdf" PDFS
rm $contig_name.depth
rm $contig_name.average-coverage
rm "windowed."$contig_name".depth"
rm $contig_name".reference_outgroup_filter.bed"
