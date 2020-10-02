# discovarif

Starting with paired reads, this pipeline will get variants such as:

- SNPs and small INDELS
- CNVs
- DELLY variants (big INDELs, translocations, duplications, inversions)

This pipeline was tested and adapted to a *Plasmodium falciparum* dataset (13 samples and 1 control) used to find drug resistance markers.

## File tree

The file tree should look like the one below.

```bash
├── out
│   ├── bambai
│   ├── fastqc
│   ├── sam
│   ├── trim
│   └── variants
│       ├── CNVs
│       │   ├── view
│       │   └── summary
│       ├── Others
│       └── SNPs-sINDELs
└── rawdata
    ├── adapters
    ├── dellysamples
    ├── Genome
    └── reads
```


Contents of this file tree are described below.

* **out**: Output directory
   * **bambai**: BAM files and their indexes (sorted and unsorted)
   * **fastqc**: FASTQC files
   * **sam**: SAM files (unsorted)
   * **trim**: Paired reads after trimming
   * **variants**: Variants
   * **CNVs**: CNV variants
     * **view**: To view CNVs in a browser (CDS, non CDS and whole regions)
     * **summary**: Filtered CNVs
   * **Others**: DELLY variants
   * **SNPs-sINDELs**: Varif variants
 * **rawdata**: Input directory
     * **adapters**: Illumina adapters used
     * **dellysamples**: DELLY tsv file for tumor and control
     * **Genome**: Must include FASTA, GFF, Core/Non core genome
     * **reads**: Sequencing reads

## Source file

A template is provided. Make sure to provide all the input files and binaries necessary for the pipeline to be launched properly.

## Pipeline

Here are quickly described the steps you can launch using the options of the `full_pipeline.sh` file.

- [Trimming/Mapping](#mapp) (*-m*)
  This option will use the raw read files and output sorted BAM files by first trimming and then mapping on the reference genome.

- [SNPs/small INDELs filtering](#varisnps) (*-s*)
  This option will use the sorted BAM files and make a VCF file with SNPs and small INDELs. The variants detected will then be filtered.

- [CNVs filtering](#varicnvs) (*-c*)

  This option will use the sorted BAM files and output per base CDS coverage in core genome for each sample. These files will then be merged and CNVs wil be filtered.

- [Other variants filtering](#varioths) (*-o*)
  This option will use the sorted BAM files and make a VCF filtered file with big INDELs, translocation, duplications and inversions.

 A thorougher description of each step is given next.

## <a name="mapp"></a>Trimming/Mapping step

This step is launched with the script `mapper-caller.sh` which can optionally run the following steps:

- [Trimming](#mapptrim) (*-q*)
- [Mapping](#mappmapp) (*-m*)
- [Variant calling](#mappvari) (*-v*)

By default, `mapper-caller.sh` will run the steps on all samples found in the reads directory, but you can subset samples by their index order using *-s*. This is useful especially when you need to parallelize the processes.

### <a name="mapptrim"></a>Trimming

Provided Illumina adapters provided in the *adapters* folder, [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) (tested with Trimmomatic v 0.39) will trim the raw reads found in the *reads* folder and output paired and unpaired read files. The paired reads will be the one used for the next steps

Additionally, the quality of raw and trimmed reads will be analysed using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (tested with FastQC v 0.10.1).

### <a name="mappmapp"></a>Mapping

The trimmed paired reads will be mapped on the reference Genome (tested with *Plasmodium falciparum* 3D7 v 46 genome) with [bwa](http://bio-bwa.sourceforge.net/bwa.shtml) (tested with bwa v 0.7.17). The output will then be processed by [samtools](http://www.htslib.org/doc/samtools.html) (tested with samtools v 1.10) in order to finally obtain sorted BAM files.

### <a name="mappvari"></a>Variant calling (each sample)

This phase has never been properly tested so it might be prone to some bugs. For each sample, [bcftools](http://samtools.github.io/bcftools/bcftools.html) will run a mpileup by reading at most 99999 reads at a position and will output a VCF file containing the DP, AD and SP fields. A variant call would then be performed using the multi-allelic model instead of the default one. Finally, only variants whose QUAL fiield is superior to 10 will be kept.

## Variant filtering

### <a name="varisnps"></a>SNPs/small INDELs

At first, bcftools (tested with bcftools v 2.26.0) will run a mpileup by reading at most 99999 reads at a position for all samples and will output a VCF file containing the DP, AD and SP fields. A variant call will then be performed using the multi-allelic model instead of the default one. Finally, only variants whose QUAL fiield is superior to 10 will be kept.

Then [varif](https://github.com/marcguery/varif) will be used to filter variants based on read depths and ALT allele frequency. A variant is considered available in a sample only if the total read depth is above 5. The ALT allele frequency of the available samples must comprise a value of at least 0.8 in one sample and a value of at most 0.2 in any other sample.

### <a name="varicnvs"></a>CNVs

After restricting the analysis on the locations of the core genome (3D7 core genome annotations provided by *[Miles A, Iqbal Z, Vauterin P, et al. Indels, structural variation, and  recombination drive genomic diversity in Plasmodium falciparum. Genome Res. 2016;26(9):1288-1299. doi:10.1101/gr.203711.115](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5052046/)*), a file containing the per base coverage in CDS regions (as annotated on the GFF file) is created for each sample using [bedtools](https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html) (tested with bedtools v 2.26.0).

Alternatively, per base coverage is generated in BED and BEDGRAPH format in order to visualise CNVs on a genome browser.

For each position in CDS regions, the coverage is divided by the median coverage of all the CDS regions. A comparison of the median normalised coverage of each CDS is made between samples in order to decide to keep it in the filtered output. 

If the median coverage of this CDS is of 1.8 or more in at least one sample and 1.5 or less in the control sample (tested with a drug sensitive sample), the CDS is selected as a potential duplication. If the median coverage of this CDS is of 0.2 or less in at least one sample and 0.5 or more in the control sample, the CDS is selected as a potential deletion.

### <a name="varioths"></a>Other variants

To find other variants such as big INDELs, duplication, translocation and inversion, [DELLY](https://github.com/dellytools/delly) (tested with DELLY v 0.7.5) somatic model will be used. Variants of all sizes will be kept if one sample has an ALT allele frequency  of 0.5 or more while the control sample has an ALT allele frequency of 0.5 at most if the read depth is at 5 or more reads.

## Setup

We tested this pipeline using the inut described below:

| Name                                                         | Version                                                      | Usage                       |
| ------------------------------------------------------------ | ------------------------------------------------------------ | --------------------------- |
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | 0.10.1                                                       | Viewing quality of reads    |
| [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) | 0.39                                                         | Trimming reads              |
| [bwa](http://bio-bwa.sourceforge.net/bwa.shtml)              | 0.7.17                                                       | Mapping reads               |
| [samtools](http://www.htslib.org/doc/samtools.html)          | 1.10                                                         | Produce BAM files           |
| [bcftools](http://samtools.github.io/bcftools/bcftools.html) | 1.10.2                                                       | Getting SNPs/small INDELs   |
| [varif](https://github.com/marcguery/varif)                  |                                                              | Filtering SNPs/small INDELs |
| [bedtools](https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html) | 2.26.0                                                       | Filtering CNVs              |
| [DELLY](https://github.com/dellytools/delly)                 | 0.7.5                                                        | Filtering other variants    |
| 3D7 genome                                                   | 46                                                           | Reference genome            |
| Core annotation                                              | [2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5052046/) | Subseting CNVs              |
| Resistant samples                                            |                                                              | Variant discovery           |
| Sensitive sample                                             |                                                              | Control sample              |

