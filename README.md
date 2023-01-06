# discovarif

Starting with paired reads, this pipeline will filter reads and get variants such as:

- SNPs and small INDELS
- CNVs
- DELLY variants (big INDELs, translocations, duplications, inversions)

This pipeline was tested and adapted to a *Plasmodium falciparum* dataset (13 samples and 1 control) used to find drug resistance markers.

The main script is launched with the command

```bash
./bin/full-pipeline -h
```

The different options available will let you select which step of the pipeline to run, the location of the configuration file and the number of processors to use.

## Tutorial

To check if you can run the pipeline, run the command

```bash
./bin/full-pipeline -b
```

Or, if you have already modified the config file:

```bash
./bin/full-pipeline -b -k config/config.sh
```

Check the output of this command to see if you have successfully installed all the required tools. This command will also tell if the folder and files that you provided in the *config* file do indeed exist.

### File tree

To work properly, the pipeline needs raw data files all located in a dedicated folder (see **DATADIR** in the config file). The entirety of the folder will be copied to a distant server if you happen to have your files remotely located (in that case files should be located under a dedicated remote folder, **REMOTEDATADIR** from the config file).

After having run all the pipeline, the file tree should look like the one below. You are only required to provide the path to the dedicated output folder (**OUTDIR** in the *config* file), the rest will be created according to the folders mentioned in the *config* file. Again if you use a remote folder, the output already created should be under a dedicated remote folder (**REMOTEOUTDIR** from the config file).

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

### `Config` file

Once you have all you raw data files located under a dedicated folder, copy the template *config-template.sh* and create your own config file.

```bash
cp config/config-template.sh config/config.sh
```

Make sure to provide all the input files and binaries necessary for the pipeline to be launched properly. Only the *TOBEFILLED* part are required.

To set up the pipeline, you should be able to provide:

* **Genome**: Must include FASTA, GFF, Core/Non core genome, ploidy
* **Reads**: The location of the reads, the sequencing adapters to be clipped

### `Samples` file

A template is provided in the `config` folder. This file should be located in the **DATADIR** folder. For example, use the command:

```
cp config/samples-template.sh $DATADIR/samples.sh
```

Each line corresponds to a sample with a uniquely identified name and both read file names mentioned in the *pair1* and *pair2* fields (the directory where the reads are located is mentioned in the `config` file). The *group* field should contain an integer equal to 0 to characterize a control sample or superior to 0 to characterize a regular group of samples. The *keep* field is used to tell the pipeline to process samples annotated as yes and exclude those annotated as no.

After the launch of the pipeline, the remaining samples to be processed (*keep* field set to *yes*) will be saved in a file whose name is identical to the file provided in the config file, with a *.run* extension prepended to the original extension.

### Run

Now that you have your own config file and your samples file, you can run the pipeline (Filtering step) using the command

```bash
./bin/full-pipeline -k config/config.sh -q
```

For more precisions n the different step, check the next paragraph.

## Pipeline

The available commands are printed with the *-h* option:

```bash
./bin/full_pipeline.sh usage:
        b) # Launch test step.
        q) # Launch quality step.
        m) # Launch mapping step.
        s) # Launch SNP/small INDEL step.
        c) # Launch CNV step.
        o) # Launch other variants step.
        k) # Use this configuration file
        n) # Process this number of samples in parallel
        u) # Launch this number of processes per sample
        h | *) # Show help.
```

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

This option uses [GATK](https://gatk.broadinstitute.org/hc/en-us) HaplotypeCaller (tested with GATK v 4.2.0) on mapped reads to obtain a GVCF file for each sample. The ploidy used by the model must be properly set in the `config` file and is by default equal to 1. Because this step can be time consuming, the GVCF files produced by the pipeline can be transferred from the temporary folder (`_tmp/gvcf.xxxxxx`) to the GVCFDIR folder from the config file. That way, whenever the pipeline detects that a sample has already a GVCF file present in the GVCFDIR, this step will not be launched.

After the [variant calling](#mappvari)  for each sample has finished, GATK CombineGVCFs will be used to merge all the GVCF and the final variants will be extracted using GATK GenotypeGVCFs.

Then [varif](https://github.com/marcguery/varif) will be used to filter variants based on read depths and ALT allele frequency. A variant is considered available in a sample only if the total read depth is above 5. Several combinations of ALT and REF allele frequencies are used to filter variants for each group as mentioned in the samples file. For example, the ALT allele frequency of at least one sample must be superior to 0.8 while being in the meantime inferior to 0.2 in at least one other sample (`--ratio-alt 0.8 --ratio-no-alt 0.2`).
Each  group will be individually filtered by varif with the control group (when a control group is provided in the samples file with samples whose *group* field is equal to 0) or without it.

### <a name="varicnvs"></a>CNVs

After restricting the analysis on the locations of the core genome (3D7 core genome annotations provided by *[Miles A, Iqbal Z, Vauterin P, et al. Indels, structural variation, and  recombination drive genomic diversity in Plasmodium falciparum. Genome Res. 2016;26(9):1288-1299. doi:10.1101/gr.203711.115](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5052046/)*), a file containing the per base coverage in CDS regions (as annotated on the GFF file) is created for each sample using [bedtools](https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html) (tested with bedtools v 2.27.1).

Alternatively, per base coverage is generated in BED and BEDGRAPH format in order to visualise CNVs on a genome browser.

For each position in CDS regions, the coverage is divided by the median coverage of all the CDS regions. A comparison of the median normalised coverage of each CDS is made between samples in order to decide to keep it in the filtered output. 

If the median coverage of this CDS is of 1.8 or more in at least 20 % of the samples and 1.5 or less in all the control samples (tested with a drug sensitive sample), the CDS is selected as a potential duplication. If the median coverage of this CDS is of 0.2 or less in at least 20 % of the samples and 0.5 or more in all the control samples, the CDS is selected as a potential deletion.

### <a name="varioths"></a>Other variants

To find other variants such as big INDELs, duplication, translocation and inversion, [DELLY](https://github.com/dellytools/delly) (tested with DELLY v 0.8.7) somatic model will be used. Variants of all sizes will be kept if one sample has an ALT allele frequency  of 0.5 or more while the control sample has an ALT allele frequency of 0.5 at most if the read depth is at 5 or more reads.

## Other options

### Multi-threading

You can choose to parallelize the processes to speed up the time required to obtain results. There are two option controlling this: -*n* for choosing the number of samples to run in parallel and -*u* for allocating a number of threads for each sample.
For example with 8 samples to process, -*n* 8 -*u* 2 will process at the same time all the 8 samples with 2 threads for each of them. 

There are optimal combinations of the -*n* and -*u* at each step to minimize the processing time; with N being equal to `min(number of samples, number of available CPUs)`:

- Filtering step:  ```-n N -u 4``` 
- Mapping: ```-n N -u 4``` 
- All variant filtering steps: ```-n N -u 1``` 

Note that some third-party tools like GATK may require a large amount of RAM which might slow down or even abort a job requesting too much processes to be launched at the same time. In that case you should reduce the number of requested processors or request more memory.

When a step of pipeline is not split by sample (like the 'Other variants' step) the number of allocated processors is the result of `(value of -n * value of -u)`. In the example above (-*n* 8 -*u* 2), this would result in 16 processors used for these steps.

### Downloading files locating in a SSH server

When raw data and output files are located on a distant server, the variables *REMOTEDATADIR*, *REMOTEADDRESS* and *REMOTEOUTDIR* must be filled. These values correspond to the location of the files on the distant server. *DATADIR* and *OUTDIR* point to the location of the files of the client. **The client must have the ability to connect to the server via a SSH key to download the files**.

The steps of the pipeline will be launched after all the files that exist in the distant server but not locally are downloaded. File modification time is not checked in the process so local files have to be removed to be replaced by those on the distant server. **The only exception to this is the `Samples` file which is always replaced by the remotely located one**.

This automatic downloading can be useful especially when the pipeline is launched via a computing grid (e.g., SGE or slurm) which requires files to be located in a node-specific folder. In that case, the *DATADIR* and *OUTDIR* should be pointing to node-specific folders, while *REMOTEDATADIR* and *REMOTEOUTDIR* should be pointing to data-specific folders.

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
| [varif](https://github.com/marcguery/varif)                  | 0.3.2                                                        | Filtering SNPs/small INDELs |
| [bedtools](https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html) | 2.26.0                                                       | Filtering CNVs              |
| [bedGraphToBigWig](https://github.com/ENCODE-DCC/kentUtils)  | 4                                                            | Compressing bed files       |
| [DELLY](https://github.com/dellytools/delly)                 | 1.1.6                                                        | Filtering other variants    |
| 3D7 genome                                                   | 46                                                           | Reference genome            |
| Core annotation                                              | [2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5052046/) | Subseting CNVs              |
| Resistant samples                                            |                                                              | Variant discovery           |
| Sensitive sample                                             |                                                              | Control sample              |

