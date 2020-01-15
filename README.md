### A pipeline written using Snakemake

*Aim: to get the value of present-day human contamination in human DNA sample*


The various tasks done by this pipeline are:

- Downloading the WGS (whole genome sequence) data files (for instance, BAM files) from the respective website

- calculating depth of coverage (doc) of data files downloaded

- downsampling the reads in order to control the contamination fraction and the doc desired

- merging endogenous and contaminants DNA followed by indexing files

- filtering data using software *angsd*
- calculating contamination rate and sequencing error using software *contaminationX*
- Plots generated using a small Python script (see **plotting.py** present in the current repository)


### Installation: ###

In order to use this pipeline, one needs to install

- Snakemake
- software *angsd* (>=0.922) 
  http://www.popgen.dk/angsd/index.php/ANGSD
- sofware *ContaminationX* 
  https://github.com/sapfo/contaminationX
- samtools
- R
  https://cran.r-project.org
- doParallel package (a R package)
- Python


### Running the pipeline: ###

* running the snakefile 
```
snakemake --cores 8
```

* to see the steps done by pipeline 

```
snakemake --dag | dot -Tpdf > name_file.pdf
```