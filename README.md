# discovarif

Starting with paired reads, this pipeline will get variants such as:

- SNPs and small INDELS
- CNVs
- DELLY variants (big INDELs, translocations, duplications, inversions)

This pipeline was tested and adapted to a *Plasmodium falciparum* dataset (13 samples and 1 control) used to find drug resistance markers.

## File tree

After having run all the pipeline, the file tree should look like the one below.

```bash
├── out
│   ├── bambai
│   ├── fastqc
│   ├── trim
│   └── variants
│       ├── CNVs
│       │   ├── view
│       │   └── summary
│       ├── Others
│       └── SNPs-sINDELs
```


Contents of this file tree are described below.

* **out**: Output directory
   * **bambai**: BAM files and their indexes (sorted and unsorted)
   * **fastqc**: FASTQC files
   * **trim**: Paired reads after trimming
   * **variants**: Variants
   * **CNVs**: CNV variants (CDS, non CDS and whole regions)
     * **view**: To view CNVs in a browser
     * **summary**: Filtered CNVs
   * **Others**: DELLY variants
   * **SNPs-sINDELs**: GATK variants

## `Config` file

A template is provided in the `config` folder. Make sure to provide all the input files and binaries necessary for the pipeline to be launched properly.

To set up the pipeline, you should be able to provide:

* **Genome**: Must include FASTA, GFF, Core/Non core genome, ploidy
* **Reads**: The location of the reads, the sequencing adapters to be clipped

## `Samples` file

A template is provided in the `config` folder. 

Each line corresponds to a sample with a uniquely identified name and both read file names mentionned in the *pair1* and *pair2* fields (the directory where the reads are located is mentionned in the `config` file). The *group* field should contain samples annotated as control or tumor. The *keep* field is used by the variant filtering step to extract only samples annotated as yes and exclude those annotated as no.

## Pipeline

Here are quickly described the steps you can launch using the options of the `full_pipeline.sh` file.

- [Filtering](#filt)  (-*q*)
  This option will remove the low quality sequences and sequencing adapters from the reads and provide a quality report.

- [Mapping](#mapp) (*-m*)
  This option will use the trimmed reads from the previous steps or if the reads provided in the `samples` file and output sorted BAM files with duplicates removed.

- [SNPs/small INDELs filtering](#varisnps) (*-s*)
  This option will use the sorted BAM files and make a VCF file with SNPs and small INDELs. The variants detected will then be filtered.

- [CNVs filtering](#varicnvs) (*-c*)

  This option will use the sorted BAM files and output per base CDS coverage in core genome for each sample. These files will then be merged and CNVs wil be filtered.

- [Other variants filtering](#varioths) (*-o*)
  This option will use the sorted BAM files and make a VCF filtered file with big INDELs, translocation, duplications and inversions.

 A thorougher description of each step is given next.

## <a name="filt"></a>Filtering step

This step is launched with the script `mapper-caller.sh` which can optionally run the trimming (*-q*) step.

With Illumina adapters provided in the `config` file, [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) (tested with Trimmomatic v 0.39) will trim the raw reads found in the *reads* folder and output paired and unpaired read files. The paired reads will be the one used for the next steps

Additionally, the quality of raw and trimmed reads will be analysed using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (tested with FastQC v 0.10.1).

## <a name="mapp"></a>Mapping step

This step is launched with the script `mapper-caller.sh` which can optionally run the mapping (-*m*) step.

The trimmed paired reads will be mapped on the reference genome (tested with *Plasmodium falciparum* 3D7 v 46 genome) with [bwa](http://bio-bwa.sourceforge.net/bwa.shtml) (tested with bwa v 0.7.17). The output will then be processed by [samtools](http://www.htslib.org/doc/samtools.html) (tested with samtools v 1.10) in order to obtain sorted BAM files. Finally, duplicated reads are removed using [picard](https://broadinstitute.github.io/picard/) (tested with picard v 2.18.25).

By default, `mapper-caller.sh` will run the steps on all samples found in the reads directory, but you can subset samples by their index order using *-s*. This is useful especially when you need to parallelize the processes. 

## Variant filtering

### <a name="varisnps"></a>SNPs/small INDELs

This step is first launched with the script `mapper-caller.sh` which can optionally run the variant calling (-*v*) step on each sample:

This option uses [GATK](https://gatk.broadinstitute.org/hc/en-us) HaplotypeCaller (tested with GATK v 4.2.0) on mapped reads to obtain a GVCF file for each sample. The ploidy used by the model must be properly set in the `config` file and is by default equal to 1.

After the [variant calling](#mappvari)  for each sample has finished, GATK CombineGVCFs will be used to merge all the GVCF and the final variants will be extracted using GATK GenotypeGVCFs.

Then [varif](https://github.com/marcguery/varif) will be used to filter variants based on read depths and ALT allele frequency. A variant is considered available in a sample only if the total read depth is above 5. The ALT allele frequency of the available samples must comprise a value of at least 0.8 in one sample and a value of at most 0.2 in any other sample.

### <a name="varicnvs"></a>CNVs

After restricting the analysis on the locations of the core genome (3D7 core genome annotations provided by *[Miles A, Iqbal Z, Vauterin P, et al. Indels, structural variation, and  recombination drive genomic diversity in Plasmodium falciparum. Genome Res. 2016;26(9):1288-1299. doi:10.1101/gr.203711.115](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5052046/)*), a file containing the per base coverage in CDS regions (as annotated on the GFF file) is created for each sample using [bedtools](https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html) (tested with bedtools v 2.27.1).

Alternatively, per base coverage is generated in BED and BEDGRAPH format in order to visualise CNVs on a genome browser.

For each position in CDS regions, the coverage is divided by the median coverage of all the CDS regions. A comparison of the median normalised coverage of each CDS is made between samples in order to decide to keep it in the filtered output. 

If the median coverage of this CDS is of 1.8 or more in at least 20 % of the samples and 1.5 or less in all the control samples (tested with a drug sensitive sample), the CDS is selected as a potential duplication. If the median coverage of this CDS is of 0.2 or less in at least 20 % of the samples and 0.5 or more in all the control samples, the CDS is selected as a potential deletion.

### <a name="varioths"></a>Other variants

To find other variants such as big INDELs, duplication, translocation and inversion, [DELLY](https://github.com/dellytools/delly) (tested with DELLY v 0.8.7) somatic model will be used. Variants of all sizes will be kept if one sample has an ALT allele frequency  of 0.5 or more while the control sample has an ALT allele frequency of 0.5 at most if the read depth is at 5 or more reads.

## Other options

### Multithreading

You can choose to parallelize the processes to speed up the time required to obtain results. There are two option controlling this: -*t* for allocating a total number of threads and -*u* for allocating a number of threads for each sample.
For example with 8 samples to process, -*t* 16 -*u* 2 will process at the same time all the 8 samples with 2 threads for each of them.
With 8 samples again, -*t* 16 -*u* 4 will process 4 samples at the same time with 4 threads for each of them.

There are specific combinations of the -*t* and -*u* at each step to minimize the processing time (with N being the number of samples):

- Filtering step:  ```-t N*4 -u 4``` 
- Mapping: ```-t N*4 -u 4``` 
- All variant filtering steps: ```-t N -u 1``` 

### Downloading files locating in a SSH server

When raw data and output files are located on a distant server, the variable *REMOTEDATADIR*, *REMOTEADDRESS* and *REMOTEOUTDIR* must be filled. These values correspond to the location of the files on the distant server. *DATADIR* and *OUTDIR* point to the location of the files of the client. The client must have the ability to connect to the server via a SSH key to download the files.

The steps of the pipeline will be launched after all the files that exist in the distant server but not locally are downloaded. File modification time is not checked in the process so local files have to be removed to be replaced by those on the distant server.

## Setup

We tested this pipeline using the programs/inputs described below:

| Name                                                         | Version                                                      | Usage                       |
| ------------------------------------------------------------ | ------------------------------------------------------------ | --------------------------- |
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | 0.10.1                                                       | Viewing quality of reads    |
| [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) | 0.39                                                         | Trimming reads              |
| [bwa](http://bio-bwa.sourceforge.net/bwa.shtml)              | 0.7.17                                                       | Mapping reads               |
| [picard](https://broadinstitute.github.io/picard/)           | 2.18.25                                                      | Removing duplicated reads   |
| [samtools](http://www.htslib.org/doc/samtools.html)          | 1.10                                                         | Produce BAM files           |
| [bcftools](https://samtools.github.io/bcftools/bcftools.html) | 1.10.2                                                       | Extract DELLY variant files |
| [GATK](https://gatk.broadinstitute.org/hc/en-us)             | 4.2.0.0                                                      | Getting SNPs/small INDELs   |
| [varif](https://github.com/marcguery/varif)                  | 0.1.1                                                        | Filtering SNPs/small INDELs |
| [bedtools](https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html) | 2.26.0                                                       | Filtering CNVs              |
| [bedGraphToBigWig](https://github.com/ENCODE-DCC/kentUtils)  | 4                                                            | Compressing bed files       |
| [DELLY](https://github.com/dellytools/delly)                 | 0.8.7                                                        | Filtering other variants    |
| 3D7 genome                                                   | 46                                                           | Reference genome            |
| Core annotation                                              | [2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5052046/) | Subseting CNVs              |
| Resistant samples                                            |                                                              | Variant discovery           |
| Sensitive sample                                             |                                                              | Control sample              |

