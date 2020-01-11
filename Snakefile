all_conta_pop = {"CEU" : ["NA12812", "NA12829", "NA12413", "NA12716"], "JPT" : ["xyz", "abc"], "CHB" : ["fgh", "mno"]} # all male samples from 1000 genomes project
# one can consider NGS data from other panels too like HGDP, SGDP etc.

def get_bam_to_merge(wildcards):
    doc_endogenous = float(wildcards.doc)*(1 - float(wildcards.conta_frac))
    doc_contaminant = (1/float(wildcards.nb_contaminant))*float(wildcards.doc)*(float(wildcards.conta_frac))
    file_endo = [f"downsampled_bam/doc_{doc_endogenous}/{wildcards.sample_endo}.bam"]
    conta_pop = all_conta_pop[wildcards.population]
    all_file_contaminant = [f"downsampled_bam/doc_{doc_contaminant}/{conta_pop[i]}.bam" for i in range(int(wildcards.nb_contaminant))]
    all_files = all_file_contaminant + file_endo
    return(all_files)

DoC = [1.00] # depth of coverage
conta = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90]
#conta = [0.01, 0.10, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99] # different values of contamination fraction
config["conta"] = conta

rule all :
     input:
        expand("plots_{doc}/{sample_endo}_{population}_{nb_contaminant}.eps",
        doc = DoC, cont = conta,
        sample_endo = ["NA19213"], # endogenous DNA is a Yoruba individual from 1000 Genomes project
        population = ["CEU"], # population id
        nb_contaminant = ["1"]  # number of contaminants
        )

rule download_sample_info:
    output:
        "20130502.phase3.low_coverage.alignment.index"
    shell:
        "wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20130502.phase3.low_coverage.alignment.index"


rule download_data:  # rule to download BAM files
    input:
        "20130502.phase3.low_coverage.alignment.index"
    output:
        "initial_bam/{sample}.bam"
    shell:
        "url=$(cut -f1 {input} | grep '\.mapped' | grep {wildcards.sample});"
        "samtools view -b ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/$url X > {output}"


rule compute_doc:    # rule to calculate the depth of coverage of bam files downloaded
    input:
        "initial_bam/{sample}.bam"
    output:
        "{sample}.doc"
    shell:
        "samtools depth -a {input} |  awk '{{sum+=$3}} END {{print sum/NR}}' > {output} "


rule downsampling_reads:  # rule to downsample reads in order to control the contamination fraction and the depth of coverage desired
    input:
        bamfile = "initial_bam/{sample}.bam",
        docfile = "{sample}.doc"
    output:
        "downsampled_bam/doc_{final_doc}/{sample}.bam"
    shell:
        "doc=$(cat {input.docfile});frac=$(echo 'scale = 5;' {wildcards.final_doc} / $doc | bc);samtools view -s $frac -b {input.bamfile} > {output}"


rule merge:  # rule for merging endogenous DNA with contaminants DNA
    input:
        get_bam_to_merge
    output:
        "mixed_bam/doc_{doc}/conta_{conta_frac}/{sample_endo}_{population}_{nb_contaminant}.bam"
    shell:
        "samtools merge {output} {input}"


rule indexing:
     input:
         "{sample}.bam"
     output:
         "{sample}.bam.bai"
     shell:
         "samtools index {input}"


rule counting_no_of_bases_part1:  # rule for filtering data using angsd software, part 1
     input:
         bam_file = "{bamfile_name}.bam",
         bam_bai_file = "{bamfile_name}.bam.bai"
     output:
         "{bamfile_name}.icnts.gz"
     shell:
         "~/bin/angsd -i {input.bam_file} -r X: -doCounts 1 -iCounts 1 -minMapQ 30 -minQ 20 -out {wildcards.bamfile_name}"
         # here, I have specified the path (~/bin) of executable file angsd (a part of angsd software)


rule counting_no_of_bases_part2:   # rule for filtering data using angsd software, part 2
     input:
         "{bamfile_name}.icnts.gz"
     output:
         "{bamfile_name}.count"
     shell:
         "~/bin/contamination -b 5000000 -c 154900000 -k 1 -m 0.05 -d 3 -e 20 -h ~/bin/contaminationX/HapMapFreqs/HapMapJPT.gz -a {input} > {output}"  # a catch here! "{input}" me!
         # here, I have specified the path (~/bin) of executable file contamination (a part of angsd software)


rule getting_contamination_and_seq_error: # rule to calculate contamination rate using contaminationX software
     input:
         count = "{bamfile_name}.count"
     output:
         "{bamfile_name}_result.txt"
     shell:
         "Rscript ~/bin/contaminationX/bin/ContaEstBoth.R counts={input.count} freqs=~/bin/contaminationX/HapMapFreqs/HapMapJPT.gz maxsites=1000 nthr=8 outfile={output} oneCns=1"
          # here, I have specified the path (~/bin) of contaminationX software


rule gather_results_files:
    input:
        expand("mixed_bam/doc_{{doc}}/conta_{conta}/{{sample_endo}}_{{population}}_{{nb_contaminant}}_result.txt", conta = conta)
    output:
        "merged_results_files/doc_{doc}/{sample_endo}_{population}_{nb_contaminant}_merged_final.txt"
    shell:
        "cat {input} > {output}"


rule getting_figures:     # rule to get the plots using a Python script
     input:
         "merged_results_files/doc_{doc}/{sample_endo}_{population}_{nb_contaminant}_merged_final.txt"
     output:
         "plots_{doc}/{sample_endo}_{population}_{nb_contaminant}.eps" # to get the plots
     script:
         "plotting.py"
