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
│       │   ├── bed
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
     * **bed**: To view CNVs in a browser (CDS, non CDS and whole regions)
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
  This option will use the sorted BAM files and make a VCF. The variants detected will then be filtered.

- [CNVs filtering](#varicnvs) (*-c*)

  This option will use the sorted BAM files and output per-base CDS coverage in core genome for each sample using [bedtools](https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html). These files will then be merged CNVs wil be filtered by `CNV-caller.R` script.

- *-o*
  This option will use the sorted BAM files and output [DELLY](https://github.com/dellytools/delly) variants in a VCF format.

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

Then [varif](https://github.com/marcguery/varif) will be used to filter variants based on read depths and $ALT\over(ALT+REF)$ ratios. A variant is considered available in a sample only if the total read depth is above 5. The $ALT\over(ALT+REF)$ ratios of the available samples must comprise a value of at least 0.8 in one sample and a value of at most 0.2 in any other sample.

### <a name="varicnvs"></a>CNVs

