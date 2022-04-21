# genomic_enrichment

### Dependencies: 
#### - python3: https://www.python.org/downloads/
#### - pybedtools: https://daler.github.io/pybedtools/main.html
#### - BEDTools: https://github.com/arq5x/bedtools2
#### - A C/C++ compiler

### Installation Steps:

#### 1. On the "main" branch, go to "Code" and select "Download ZIP" from the dropdown menu and place the zip file in a directory where it can be extracted.

#### 2. After succesfully extracting from the zip, enter the "genomic_enrichment-main" directory that has been newly created, the directory should have the python script "calculate_enrichment.py".

#### 3. On the local machine command line interface, while having all the required dependencies installed, execute `python calculate_enrichment.py --help` and the options help menu should be printed to the screen.

#### 4. Please create an additional directory in the same directory level as the python script named "genomeFASTA" and place your genome Fasta files in there (genomeFASTA/hg.19).

### Usage:
```
usage: calculate_enrichment.py [-h] [-i ITERS] [-s {hg19,hg38,mm10,dm3}]
                               [-b BLACKLIST] [--GC_blacklist GC_BLACKLIST]
                               [-n NUM_THREADS]
                               [--print_counts_to PRINT_COUNTS_TO]
                               [--elem_wise] [--by_hap_block] [--GC_option]
                               [--GC_max GC_MAX] [--GC_min GC_MIN]
                               [--GC_margin GC_MARGIN]
                               [--GC_bp_resolution GC_BP_RESOLUTION]
                               region_file_1 region_file_2

Calculate enrichment between bed files.

positional arguments:
  region_file_1         bed file 1 (shuffled)
  region_file_2         bed file 2 (not shuffled)

optional arguments:
  -h, --help            show this help message and exit
  -i ITERS, --iters ITERS
                        number of simulation iterations; default=100
  -s {hg19,hg38,mm10,dm3}, --species {hg19,hg38,mm10,dm3}
                        species and assembly; default=hg19
  -b BLACKLIST, --blacklist BLACKLIST
                        custom blacklist file; default=None
  --GC_blacklist GC_BLACKLIST
                        custom blacklist file for GC_option (content of file
                        is up to user); default=None
  -n NUM_THREADS, --num_threads NUM_THREADS
                        number of threads; default=SLURM_CPUS_PER_TASK or 1
  --print_counts_to PRINT_COUNTS_TO
                        print expected counts to file
  --elem_wise           perform element-wise overlaps; default=False
  --by_hap_block        perform haplotype-block overlaps; default=False
  --GC_option           perform shuffling with regions of similar GC content;
                        default=False
  --GC_max GC_MAX       custom max GC percent threshold (integer, ex: 80) for
                        GC_option, must be used with --GC_min, default is
                        --GC_margin default settings
  --GC_min GC_MIN       custom min GC percent threshold (integer, ex: 20) for
                        GC option, must be used with --GC_max, default is
                        --GC_margin default settings
  --GC_margin GC_MARGIN
                        adjust GC content allowed margin in GC_option;
                        default=0.1(10% GC content error margin)
  --GC_bp_resolution GC_BP_RESOLUTION
                        adjust GC content bp resolution in GC_option;
                        default=100(bp) 
```
