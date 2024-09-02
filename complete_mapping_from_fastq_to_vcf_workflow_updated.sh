##################################### MAPPING FROM FASTQ to VCF WORKFLOW ##########################################################################
##########################################################################################################################################


###########################################################################################################################################
#STEP 0.5: Get the fastq files from ENA
#get SRR number for the accessions of interest
https://www.ebi.ac.uk/ena/browser/view/PRJNA273563
	#PRJNA273563: this number has to be known from people (publications) who submitted the data 
	#note that accession ID we see in 1001genome database is called "sample_alias" in ENA
	#SRR number is run_accession

#install sra toolkit to use "prefetch" and "fastq-dump" to download the fastq files using SRR number 
    #add programme to your bash profile
export PATH=/path/to/your/directory/sratoolkit.3.0.5-centos_linux64/bin:$PATH

#prefetch first to get sra files # this can even be done without submitting jobs to CHEOPS (336 samples took a night)
prefetch $(</path/to/your/directory/SRA_list.txt) -O /path/to/your/directory/srafiles

#then "convert" the fastq files from SRA files with fastq-dump
#BIG NOTE: MUST run in the directory where your SRA folders are, or it won't work 
#this can take forever so better submit ut to cheops. See the fastq_from_sra.sh script
fastq-dump --split-3 --gzip --outdir /path/to/your/directory/fastqfiles $(</path/to/your/directory/SRA_list.txt) 

##########################################################################################################################################
#STEP 1: Use fastp to filter low quality signal, remove polyG and adapters 
#(TRIMMOMATIC needs specification but fastp can automatically detect the adpapters)
    #non paired end works the same, just need to specify --in1 and --out1 then

mkdir filteredReads
mkdir mappedReads
    #try to be organised with raw reads, filtered reads, and later mapped reads like this, or you might go insane in the end

#now we are at the raw read folder
for R1 in SRR*_1.fastq.gz; do
    R2=$(echo $R1 | sed 's/_1.fastq.gz/_2.fastq.gz/')
    base=$(basename "$R1" _1.fastq.gz)
    fastp --in1 "$R1" --in2 "$R2" --out1 "$base.trimmed.R1.fastq.gz" --out2 "$base.trimmed.R2.fastq.gz" -l 50 -h "$base.html" &> "$base.log"
done

#move your stuff to new folder
mv *trimmed.*.fastq.gz /path/to/your/directory/filteredReads

#it will iterate over all paired-end files, run fastp on each pair of files, 
#and write the trimmed output files with the base name of the input files followed by ".trimmed.R1.fastq.gz" and ".trimmed.R2.fastq.gz". 
#The sed command is used to replace the _1.fastq.gz with _2.fastq.gz to get the corresponding R2 file. 
#The -l 50 flag sets the minimum length of reads to 50bp, and -h writes the HTML report of quality control of fasta files. 
#The log file for each pair of files will be named using the base name of the input files followed by ".log".

#parameters to change if you wish, I used the default options
#-q : Quality threshold per base required. 
    #Default: 15, which means that a Phred quality score of at least 15 is required
#-u : Percent of bases allowed to be below the quality threshold to keep the read (0~100). 
    #Default 40 means 40% bases can fail the quality threshold. If more bases fail, the read is removed.


##########################################################################################################################################
#STEP 2: map reads
#you need to have bwa, samtools, and bcftools for the next steps 
#index reference genome (let say it's named TAIR10.fasta)
bwa index TAIR10.fasta

#make a variable pointing to the reference
REF=/path/to/your/directory/TAIR10.fasta
TRIM_DIR=/path/to/your/directory/filteredReads
OUT_DIR=/path/to/your/directory/mapppedReads


find "${TRIM_DIR}" -name "*.trimmed.R1.fastq.gz" | while read -r R1; do
    R2="${R1%.trimmed.R1.fastq.gz}.trimmed.R2.fastq.gz"
    SAMPLE=$(basename "${R1%.trimmed.R1.fastq.gz}")
    bwa mem -M -t 16 "${REF}" "${R1}" "${R2}" > "${OUT_DIR}/${SAMPLE}.sam"
done

## if needed: looking for multiple patterns with -o flag (-o stands for "or")
find "${TRIM_DIR}" \( -name "SRR19458*.trimmed.R1.fastq.gz" -o -name "SRR1946*.trimmed.R1.fastq.gz" \) | while read -r R1; do
    R2="${R1%.trimmed.R1.fastq.gz}.trimmed.R2.fastq.gz"
    SAMPLE=$(basename "${R1%.trimmed.R1.fastq.gz}")
    bwa mem -M -t 16 "${REF}" "${R1}" "${R2}" > "${OUT_DIR}/${SAMPLE}.sam"
done


#loop through all the _1.trimmed.fastq.gz files in the TRIM_DIR directory. 
#For each _1.trimmed.fastq.gz file, you extract the corresponding _2.trimmed.fastq.gz file, 
#as well as the sample name (which is the base name of the _1.trimmed.fastq.gz file). 
#You then run bwa mem on the pair of trimmed fastq files and 
#output the SAM file to the OUT_DIR directory with the sample name as the file name.

##########################################################################################################################################
#STEP 3: Convert SAM to BAM to save space and sort the BAM files
#to do it in one go
for file in *_1.trimmed.fastq.gz; do
    base=$(basename ${file%_1.trimmed.fastq.gz})
    bwa mem -M -t 16 $REF \
        ${base}_1.trimmed.fastq.gz \
        ${base}_2.trimmed.fastq.gz \
        > ${base}.sam
    samtools view -@ 12 -bS ${base}.sam > ${base}.bam
done


## IMPORTANT: before you proceed, you must SORT THE BAM FILE
## Most downstream analyses require this in order to function properly

### Sorting is done in two steps
## Loop through each .sam file to convert it to bam and sort it

cd mappedReads
for file in *.sam; do
    base=$(basename ${file%.sam})
    samtools view -@ 12 -bS ${base}.sam > ${base}.bam
done

for bamfile in *.bam
do
    sorted_bam_base=$(basename ${bamfile%.bam})_sort.bam
    
    samtools sort "$bamfile" -o "$sorted_bam_base"
done


#For each file, we generate the output filename for the sorted BAM file by removing the .sam extension from the input filename and adding _sort.bam. 
#Finally, we use samtools sort to sort the input SAM file and save the output as the sorted BAM file with the generated filename.

##########################################################################################################################################
#STEP 3.5: Optional depending on the data - Mark duplicates
## DO keep the order of the four samtools steps 
for file in *_sort.bam; do
    bamquery_sort=$(basename "$file" _sort.bam)_querysort.bam
    fixmate_output=$(basename "$file" _sort.bam)_fixmate.bam
    fixmate_sort=$(basename "$fixmate_output" _fixmate.bam)_fixmate_sort.bam
    output=$(basename "$file" _sort.bam)_sort_rmd.bam
    samtools sort -n -o "$bamquery_sort" "$file"
    samtools fixmate -m -@ 12 -O bam "$bamquery_sort" "$fixmate_output"
    samtools sort "$fixmate_output" -o "$fixmate_sort"
    samtools markdup -r -s "$fixmate_sort" "$output"
    rm "$fixmate_output" "$fixmate_sort"
done

##########################################################################################################################################
#STEP 3.5: IMPORTANT - check quality of the mapping
for file in *_sort.bam; do
    sample_name=$(basename "$file" _sort.bam)
    output="${sample_name}.stats"
    samtools stats -@ 12 "$file" > "$output"
done
    #this will create .stats file, which can be visualised by multiQC as demonstrated below
    #you should give the output .stats extension, otherwise, it will only generate some .out files 

### easily visualise with multiQC
### remember to export the PATH to anaconda (for me I installed it by conda install multiqc)
multiqc /path/to/directory/where_the_stats_files_are/ /also/path/to/directory/where_the_fastq_log_files_are/
    #multiqc will search for all stats and log files in the given directories to generate a report for all samples
    #you can give as many directories as you like

    #in case you might need this note for installing multiqc with conda on the cluster
    #in case it didn't work with the first try, do this first
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    #then try to install multiqc again
    conda install multiqc


##########################################################################################################################################
#STEP 4: Variant calling using bcftools

#first, indexing the reference again! 
#This is not redundant, the indexing by bwa is for mapping, and this actually has to be done with samtools. 

samtools faidx TAIR10.fasta
    #this create fa.fai file
    #for clarity sake, note that the columns are: chromsome name, chromosome length, offset of the first base, fasta line length, fasta line length + 1.

#first use the samtools mpileup tool to pileup our BAM files. 
#then take all reads at a given position and call variants from the reads covering that position. 
#After performing the pileup, we than pass the output to bcftools call which will actually call variants. 

#most people use GATK, which I never used but it presumably has quality control during variant calling
#mpileup does not, so I use very "universal" threshold for quality control here -q 20 and -Q 30 #see below for details

#but first, to have the vcf file with the sample names I want, I need to prepare the txt file of the sample names
#the sample name of the vcf will be the SRR* part of the SRR*_sort.bam files, so:
sample_names_file="samplenames.txt" #create an output file
> "$sample_names_file" 
    #the ">" is to make sure the file is initially empty, if you're sure it's empty, just skip that

#putting stuff to the samplenames.txt file
for file in *_sort.bam
do
    sample_name=$(basename "$file" _sort.bam)
    echo "$sample_name" >> "$sample_names_file"
done 

#Finally, putting it to variant calling step
bcftools mpileup -a AD,DP,SP -q 20 -Q 30 -Ou -f $REF \
./align/*_sort.bam | bcftools call -v -f GQ,GP \
-m -Oz -o ./outputname.vcf.gz

#basically the vcf is to be created first then I use reheader to rename all the samples with names I want
bcftools reheader --samples samplenames.txt ./outputname.vcf.gz -o ./outputname_rename.vcf.gz

#And this is the vcf file that are ready to proceed with exciting analyses!

#Flags of bcftools mpileup:


#-v only output variant sites
#-a - Annotate the vcf - here we add allelic depth (AD), genotype depth (DP) and strand bias (SP).
#-O - the output type. Here it is u which means we do not compress the output.
#-f - specify the reference genome to call variants against.
#-q - minimum mapping quality for an alignment to be used
#-Q - minimum quality for a base to be considered

#Flags of bcftools call:
#-f - format fields for the vcf - here they are genotype quality (GQ) and genotype probability (GP).
#-v - output variant sites only - i.e. ignore non-variant parts of the reads #which I did not use
#-m- use bcftools multiallelic caller
#-O- specify the output type, here it is z - i.e. gzipped (compressed) vcf
#-o output path

### Further options

#--min-mapping-quality: Specifies the minimum mapping quality threshold for reads to be considered for variant calling.
#--min-base-quality: Specifies the minimum base quality threshold for bases to be considered for variant calling.
#--min-alternate-fraction: Specifies the minimum alternate allele fraction for a site to be considered variant.
#--genotype-qualities: Includes genotype quality scores in the output VCF file.
#--use-best-n-alleles: Specifies the maximum number of alleles to output per site.
#--report-monomorphic: Includes monomorphic sites in the output VCF file.
#--allele-balance-priors-off: Disables allele balance priors.
#--allele-balance-priors-strict-require-AF: Requires the alternate allele frequency to be above 0 for allele balance priors.
#-F: Specifies the maximum fraction of sites within a population that can be missing data for a variant to be reported.
#-C: Specifies the minimum per-sample coverage required to consider a variant.
#--pooled-continuous: Enables calling in pooled samples with continuous ploidy.
#--skip-coverage: Skips variant calling for regions with a coverage higher than the specified threshold.
#--skip-read-groups: Skips outputting per-read-group information.

##########################################################################################################################################
#################################FILTERING VCF#######################################################################

#Depth: always include a minimum depth filter and ideally also a maximum depth one too. 
    #Minimum depth cutoffs will remove false positive calls and will ensure higher quality calls as well. 
    #A maximum cut off is important because regions with very, very high read depths are likely repetitive ones mapping to multiple parts of the genome.
#Quality Genotype quality is also an important filter - 
    #essentially you should not trust any genotype with a Phred score below 20 which suggests a less than 99% accuracy.
#Minor allele frequency MAF can cause big problems with SNP calls - 
    #and also inflate statistical estimates downstream. 
        #Ideally you want an idea of the distribution of your allelic frequencies but 0.05 to 0.10 is a reasonable cut-off. 
        #some analyses, particularly demographic inference can be biased by MAF thresholds.
#Missing data: How much missing data are you willing to tolerate? 
    #It will depend on the study but typically any site with >25% missing data should be dropped.

#First you need to randomly sample your vcf for the statistics, big vcf would take forever
#randomness can be "reliably" achieved by 
bcftools view outputname_rename.vcf.gz | vcfrandomsample -r 0.012 > outputname_rename_subset.vcf
#-r means the fraction of variants we want to retain
#we just want a fraction to see how things are in our vcf

# compress and index subset vcf
bgzip outputname_rename_subset.vcf
bcftools index outputname_rename_subset.vcf.gz


#Now let's get the statistics
SUBSET_VCF=/path/to/your/directory/outputname_rename_subset.vcf.gz
OUT=/path/to/your/directory/outputname_rename_subset

#calculate allele frequency for each variant
    #freq2 outputs the frequencies without information about the alleles
    #add max-alleles 2 to exclude sites that have more than two alleles.
vcftools --gzvcf $SUBSET_VCF --freq2 --out $OUT --max-alleles 2

#calculate the mean depth of coverage per individual and per site respectively
vcftools --gzvcf $SUBSET_VCF --depth --out $OUT
vcftools --gzvcf $SUBSET_VCF --site-mean-depth --out $OUT

#extract the site quality score for each site
vcftools --gzvcf $SUBSET_VCF --site-quality --out $OUT

#another important step, calculate the proportion of missing data per sample and per site respectively
vcftools --gzvcf $SUBSET_VCF --missing-indv --out $OUT
vcftools --gzvcf $SUBSET_VCF --missing-site --out $OUT

######################
### Now bring them to R
library(tidyverse)
    #first, variant quality
var_qual = read_delim("./outputname_rename_subset.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
        
a = ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light()
#Phred (qual) score of 30 represents a 1 in 1000 chance that SNP call is erroneous. 
    #so if most sites exceed this - suggesting we have a lot of high confidence calls. 

    #variant mean depth
var_depth = read_delim("./outputname_rename_subset.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
summary(var_depth$mean_depth)
b = ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() #+ xlim(0, 100)
#10x is a good rule of thumb as a minimum cutoff for read depth, 
    #although if we wanted to be conservative, we could go with 15x.

#What is more important here is that we set a good maximum depth cufoff. 
    #if some regions clearly have extremely high coverage and 
        #this likely reflects mapping/assembly errors and also paralogous or repetitive regions. 
        #We want to exclude these as they will bias our analyses. 
    #Usually a good rule of thumb is something the mean depth x 2 - so in this case we could set our maximum depth at 40x.

    #variant missingness
var_miss = read_delim("./outputname_rename_subset.lmiss", delim = "\t",col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
summary(var_miss$fmiss)
d = ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light()

#Conservative threshold is to remove all sites where over 10% of individuals are missing a genotype. 
    #NOTE here that vcftools inverts the direction of missigness, so our 10% threshold means we will tolerate 90% missingness
    #be careful when threshold is 75% etc...

    #minor allele frequency
var_freq = read_delim("./outputname_rename_subset.frq", delim = "\t", col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
# find minor allele frequency
var_freq$maf = var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
summary(var_freq$maf)
e = ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light()

#For clarity: With 16 individuals, there are 28 alleles for a given site. 
    #so f.e., MAF = 0.04 is equivalent to a variant occurring as one allele in a single individual (i.e. 28 * 0.04 = 1.12). 
    #Alternatively, an MAF of 0.1 would mean that any allele would need to occur at least twice (i.e. 28 * 0.1 = 2.8).
# as setting MAF threshold is not straightforward, it is best practice to produce one dataset with a good MAF threshold and keep another without any MAF filtering at all. 
#Common way is to set our MAF to 0.1

    #now individual statitics
#mean depth per individual
ind_depth = read_delim("./outputname_rename_subset.idepth", delim = "\t",col_names = c("ind", "nsites", "depth"), skip = 1)
f = ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light()
#misisng data
ind_miss  = read_delim("./outputname_rename_subset.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
g = ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light()


##############################################
##OKAY now we know what threhold to choose
VCF_IN=~/vcf/outputname_rename.vcf.gz
VCF_OUT=~/vcf/outputname_rename_filtered.vcf.gz

# set filters
MAF=0.05
MISS=0.95
QUAL=30
MIN_DEPTH=10
MAX_DEPTH=50

#### perform the filtering with vcftools
vcftools --gzvcf $VCF_IN \
--remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > \
$VCF_OUT


###########################################################################################################################################
########################################## Various POST PROCESSING ########################################################################

##########################################
####visualise MAF
#create MAF file
plink --vcf input.vcf.gz --freq --allow-extra-chr --out output_maf_prefix
	#allow-extra-chr to force it to proceed in case of mitochondrial chromosome f.e. 
	#visualise in R by
h = hist(freq$MAF)
plot(h, freq = FALSE)
	

##########################################
###convert vcf to plink format, bed file
plink --vcf input.vcf.gz --make-bed --out prefix_of_plink_bfiles
	#creating three output files named prefix_of_plink_bfiles


##########################################
###LD decay calculation
#with plink
plink --bfile input.vcf --r2 --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0.2 --out myld_output_prefix

#visualise LD using PopLDdecay. This takes a while but very nice
PopLDdecay-3.42/bin/PopLDdecay -InVCF input.vcf.gz -OutStat output_prefix
perl PopLDdecay-3.42/bin/Plot_OnePop.pl -inFile output_prefix.stat.gz -output output_prefix_FIGURE

##########################################
#########LD pruning and PCA
#much faster with plink than bcftools

#LD pruning with bcftools
bcftools +prune -m 0.2 -w 50000 input.vcf.gz -Oz -o output_ldprune_prefix

#LD pruning with plink, recommended, way faster than bcftools
plink --vcf input --indep-pairwise 50 10 0.1 --allow-extra-chr --out output_ldprune_prefix
	#creating 2 files named output.prune.in and output.prune.out # the first one fell below the LD threshold which should be retained, and vice versa
	#the indep-pairwise is linkage pruning, here means 50kb window, 10 is window step size, r squared threshold is set to 0.1, pruning any variables with r squared greater than 0.1
	
#perform PCA
plink --vcf input.vcf.gz --extract output_ldprune_prefix.prune.in --make-bed --pca --out output_pca_prefix

	#--make-bed is to write additional files of another type of population structure analysis - a model based approach with admixture
	# pca output is output.eigenval and output.eigenvec # also binary output output.bed for admixture analysis
	# output.bim map file of variants contained in bed file # output.fam map file of individuals contained in bed file

    #proceed with these bed, bim, fam files for gwas later

######## plot PCA in R ########
########  see draw_PCA_for_plink_output.R script #####


################ perform GWAS with gemma ################
##########################################################
### Before gwas, extract samples you want to work with
bcftools input.vcf.gz -S samplelist.csv -Oz -o output.vcf.gz
## extract chr 1 to 5 only
plink --bfile input_prefix --chr 1-5 --allow-extra-chr --make-bed --out filtered_plink_prefix
## make plink files
plink --vcf input --make-bed --out filtered_plink_prefix

#to get map and ped files (if needed)
plink --bfile prefix_output --recode --tab --out prefix_output

######GWAS######
#Run GWAS with Plink bfile --highly recommended

#### create kinship first
    #by plink
plink --bfile prefix_of_plink_bfiles --cluster --matrix --out your_kinship_output_name
    #by gemma - recommended
gemma --bfile prefix_of_plink_bfiles -p phenotype_file.txt -gk -o your_kinship_output_name

####running gwas (-n 1 is to choose phenotype 1 if data contains multiple phenotypes)
###NOTE! that the phenotype.txt file format must be pheno1, pheno2, etc 
## Make sure the order is correct since there is no sample name info here!

gemma --bfile prefix_of_plink_bfiles -p phenotype_file.txt -n 1 -k your_kinship_output_name.cXX.txt -lmm 4 -o your_gwas_output_name

###bring the file to R for manhattan plot with qqman package, manhattan function

############################## Good Luck!! ####################################
