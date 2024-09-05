# Define dataset
dataset = [
    "M1_LA", "M1_RA", "M2_LA", "M2_RA", "M3_LA", "M3_RA", "M4_LA", "M4_RA",
    "M5_LA", "M5_RA", "M6_LA", "M6_RA", "M7_LA", "M7_RA", "M8_LA", "M8_RA",
    "M9_LA", "M9_RA", "M10_LA", "M10_RA", "M11_LA", "M11_RA", "M12_LA", "M12_RA",
    "M13_LA", "M13_RA", "M14_LA", "M14_RA", "M15_LA", "M15_RA", "M16_LA", "M16_RA",
    "M17_LA", "M17_RA", "M18_LA", "M18_RA", "M19_LA", "M19_RA", "M20_LA", "M20_RA",
    "M22_LA", "M22_RA", "M23_LA", "M23_RA", "M24_LA", "M24_RA", "M25_LA", "M25_RA", 
    "Dorado_LA", "Im_A_Mets_Fan_LA", "Jytte_LA", "Kevin_Cook_LA", "San_Diego_LA", "Styles_LA",
    "Dorado_RA", "Im_A_Mets_Fan_RA", "Jytte_RA", "Kevin_Cook_RA", "San_Diego_RA", "Styles_RA"]

################################################################################
################################################################################


# Assess quality of reads with FastQC
rule unzip_fastq:
    input:
        "{folder}-data/{file}.fq.gz"
    output:
        temp("{folder}-data/{file}.fastq")
    shell:
        """
        gunzip -c {input} > {output}
        """


rule fastqc:
    input:
        "{folder}-data/{file}.fastq"
    output:
        html = "fastqc/{folder}-data/{file}.html",
        zip = "fastqc/{folder}-data/{file}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    log:
        "logs/fastqc/{folder}/{file}.log"
    threads: 1
    resources:
        mem_mb = 1024
    wrapper:
        "v3.3.3/bio/fastqc"


# MultiQC report for the raw data only
rule raw_data:
    input:
        expand("fastqc/raw-data/{sample}_R1_fastqc.zip", sample = dataset),
        expand("fastqc/raw-data/{sample}_R2_fastqc.zip", sample = dataset)
    output:
        "qc-raw-report/multiqc.html",
        directory("qc-raw-report/multiqc_data"),
    log:
        "logs/multiqc-raw.log",
    wrapper:
        "v3.3.3/bio/multiqc"

################################################################################
################################################################################

# pre-process reads with fastp
# https://www.ncbi.nlm.nih.gov/pmc/articles/pmc6129281/
# please refer to handbook (https://github.com/opengene/fastp) for details on the steps
rule fastp_pe:
    input:
        sample = ["raw-data/{sample}_R1.fq.gz", "raw-data/{sample}_R2.fq.gz"]
    output:
        trimmed = ["clean-data/{sample}_1.fq.gz", "clean-data/{sample}_2.fq.gz"],
        #failed = "qc-rejected/{sample}.failed.fastq.gz",
        html = "clean-data/{sample}.html",
        json  = "clean-data/{sample}-fastp.json"
    log:
        "logs/fastp/{sample}.log"
    params:
        # parameters for processing
        extra=" ".join([
            # average quality score of reads requirement
            "--average_qual=20",
            # at most 10% of read positions have score below 20
            "--qualified_quality_phred=20",
            "--unqualified_percent_limit=10",
            # ensure a minimal read length
            "--length_required=40",
            #  illumina truseq adapters
            "--adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            "--adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
        ])
    threads: 2
    wrapper:
        "v3.3.3/bio/fastp"


# MultiQC report for the cleaned data
rule clean_data:
    input:
        expand("fastqc/clean-data/{sample}_1_fastqc.zip", sample = dataset),
        expand("fastqc/clean-data/{sample}_2_fastqc.zip", sample = dataset),
        expand("clean-data/{sample}-fastp.json", sample = dataset)
    output:
        "qc-clean-report/multiqc.html",
        directory("qc-clean-report/multiqc_data"),
    log:
        "logs/multiqc-clean.log",
    wrapper:
        "v3.3.3/bio/multiqc"


################################################################################
################################################################################


# index the genome, but unzip first because of star (genome and annotation should be in a subfolder called data)
rule unzip_genome:
    input:
        "data/{genome}.genome.fa.gz"
    output:
        temp("data/{genome}.fa")
    shell:
        """
        gunzip -c {input} > {output}
        """


rule star_index:
    input:
        fasta = "data/{genome}.fa"
    output:
        directory("{genome}_index")
    threads: 8
    log:
        "logs/star_index_{genome}.log"
    wrapper:
        "v3.3.3/bio/star/index"

# map cleaned reads against the genome INSERT CORRECT NAME FOR GENOME
rule star_pe_multi:
    input:
        fq1 = "clean-data/{sample}_1.fq.gz",
        fq2 = "clean-data/{sample}_2.fq.gz",
        idx = "Equus_caballus.EquCab3.0.dna.toplevel_index",
    output:
        # see star manual for additional output files
        aln = "star/{sample}.bam",
        log = "star/{sample}.Log.out",
        log_final = "star/{sample}.Log.final.out",
        sj  ="star/{sample}.sj.out.tab",
    log:
        "logs/star/{sample}.log"
    threads: 8
    wrapper:
        "v3.3.3/bio/star/align"


# report on the mapping
rule map_data:
    input:
        expand("star/{sample}.Log.final.out", sample = dataset)
    output:
        "mapping-report/multiqc.html",
        directory("mapping-report/multiqc_data"),
    log:
        "logs/multiqc-mapping.log",
    wrapper:
        "v3.3.3/bio/multiqc"


################################################################################
################################################################################
#sort bam files for qualimap
rule samtools_sort:
    input:
        "star/{sample}.bam",
    output:
        "star/{sample}.sorted.bam",
    log:
        "logs/samtoolsort/{sample}.log",
    params:
        extra="-m 4G",
    threads: 8
    wrapper:
        "v3.3.5/bio/samtools/sort"

#checking quality of mapping
rule qualimap:
    input:
        # BAM aligned, splicing-aware, to reference genome
        bam="star/{sample}.sorted.bam",

    output:
        directory("qc-qualimap/{sample}"),
    log:
        "logs/qualimap/bamqc/{sample}.log",
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as {resources.mem_mb}:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=4096,
    wrapper:
        "v3.3.5/bio/qualimap/bamqc"


################################################################################
################################################################################

# unzip gff3 (annotation file, which should also be in the data-subfolder)
rule unzip_annotation:
    input:
        "data/{anno}.gff3.gz"
    output:
        temp("data/{anno}.gff3")
    shell:
        """
        gunzip -c {input} > {output}
        """

# Rule for feature counts with reverse stranded
rule feature_counts_reverse_stranded:
    input:
        samples = "star/{sample}.bam",
        annotation = "data/Equus_caballus.EquCab3.0.111.gtf"
    output:
        multiext(
            "featureCount/countReadPairs/{sample}",
            ".featureCounts",
            ".featureCounts.summary"
        )
    threads: 2
    params:
        strand = 2,  # Set to 2 for reversely stranded
        extra = " --largestOverlap -p -B -P -t exon -F GFF --countReadPairs"
    log:
        "logs/reverse_stranded/{sample}.log"
    wrapper:
        "v3.3.3/bio/subread/featurecounts"

# Rule to combine feature counts into a single expression matrix
rule feature_counts_combine_reverse_stranded:
    input:
        expand("featureCount/countReadPairs/{sample}.featureCounts", sample=dataset)
    output:
        'gene-expression-all-reverse-stranded-countReadPairs.tsv'
    shell:
        """
        tmp=$(mktemp -d)
        first_file={input[0]}
        echo "First file to process: $first_file"
        # Check if the first file exists before using tail
        if [ -f $first_file ]; then
            tail -n +2 $first_file | cut -f 1-6 > $tmp/00_annot
        else
            echo "First file $first_file does not exist. Exiting."
            exit 1
        fi

        for i in {input}; do
            echo "Checking file: $i"
            if [ -f $i ]; then
                echo "Processing file: $i"
                bsn=$(basename $i .featureCounts)
                echo "Sample name: $bsn"
                echo $bsn > $tmp/$bsn
                tail -n +3 $i | cut -f 7 >> $tmp/$bsn
            else
                echo "File $i does not exist. Skipping."
            fi
        done
        echo "Combining files..."
        paste $tmp/* > {output}
        echo "Combination complete."
        rm -rf $tmp
        """

# Rule for MultiQC report on quantification
rule count_data:
    input:
        expand("featureCount/countReadPairs/{sample}.featureCounts.summary", sample=dataset)
    output:
        "count-report/multiqc.html",
        directory("count-report/multiqc_data")
    log:
        "logs/multiqc-count.log"
    wrapper:
        "v3.3.3/bio/multiqc"

################################################################################
################################################################################
# build a single count matrix

rule feature_counts_combine:
    input:
        expand("featureCount/reverse_stranded/{sample}.featureCounts",
               sample = dataset)
    output:
        'gene-expression-all-reverse-stranded.tsv'
    shell:
        """
        tmp=$(mktemp -d)
        # Extract the gene names etc
		# (tail excludes first line [start display from 2nd line],
		#  then cut extracts the first 6 columns)
        tail +2 {input[0]} | cut -f 1-6 > $tmp/00_annot
        # for each file extract only the counts
        for i in {input} ; do
			# the sample name from the filename without the file extension
            bsn=$(basename $i .featureCounts)
			# save sample name in a temporary file
            echo $bsn > $tmp/$bsn
			# exclude first 2 rows and only extract the last column
			# append expression counts to the temporary file
            tail +3 $i | cut -f 7 >> $tmp/$bsn
        done
        # 'paste' columns together
		# the '00_annot' file fill come first, due to the alpha-numeric sotring
        paste $tmp/* > {output}
        rm -rf $tmp
        """


################################################################################
################################################################################


# Add these rules at the end of your Snakefile

# Rule for collecting BAM files
rule collect_bam:
    input:
        expand("star/{sample}.sorted.bam", sample=dataset)
    output:
        "bam_collection.txt"
    run:
        with open(output[0], "w") as f:
            f.write("\n".join(input))

# Rule for collecting Qualimap results
rule collect_qualimap:
    input:
        expand("qc-qualimap/{sample}", sample=dataset)
    output:
        "qualimap_collection.txt"
    run:
        with open(output[0], "w") as f:
            f.write("\n".join(input))

# Modify the 'all' rule to include these collections
rule all:
    input:
        "qc-raw-report/multiqc.html",
        "qc-clean-report/multiqc.html",
        "mapping-report/multiqc.html",
        "count-report/multiqc.html",
        'gene-expression-all-reverse-stranded-countReadPairs.tsv',
        "bam_collection.txt",
        "qualimap_collection.txt"
