# longHam

LongHam is an assembly pipeline designed to assemble small fungal genomes based on a combination of long and short reads. The pipeline proceeds as follows:

![alt text](https://github.com/Gabaldonlab/longHam/blob/main/pipeline.png "Pipeline summary")

## Requirements

LongHam requires the instalation of several programs. Most of them can be installed using conda.

- Trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic)
- CANU (https://github.com/marbl/canu)
- Platanus (http://platanus.bio.titech.ac.jp/)
- DBG2OLC (https://github.com/yechengxi/DBG2OLC/tree/master/compiled)
- Sparse (https://github.com/yechengxi/DBG2OLC/tree/master/compiled)
- Masurca (https://github.com/alekseyzimin/masurca)
- BWA (http://bio-bwa.sourceforge.net/)
- Samtools (http://www.htslib.org/)
- Pilon (https://github.com/broadinstitute/pilon)
- WTDBG2 (https://github.com/ruanjue/wtdbg2)
- Ragout (https://github.com/fenderglass/Ragout)

Full paths to the different programs need to be added into the first few lines of the pipeline script (longHam.py). In addition you need to add the path to the adapters file needed for trimmomatic to be used and to the masurca master config file provided with the script.

## Execute script

Once all paths are properly changed, you can execute longHam as follows:

`python3 longHam.py -s1 shortReads1 -s2 shortReads2 -n longReads -o outFolder -t 48 -g 12000000`

Where 
- -s1 and -s2 are the files that contain raw short read data.
- -n is the file that contains long reads. These reads should have been previously cleaned
- -o Name for the outfolder where results will be stored
- -t indicates the number of threads
- -g indicates the genome size in bp

Other options

- -c The program will only use a subset of long reads to run most programs. By default this value is set to the subset that gives a 30X genome coverage
- --reads_type by default the program assumes reads are from Oxford Nanopore, if they are pacbio you should use this option
- -r In the last step the user can choose to provide a reference genome independently assembled. This will only be used by Ragout to provide a more complete assembly.

<i> All paths included should be full paths, relative paths will not work </i>

## Outputs

The main results of longHam.py will be stored in the assemblies folder within the folder you provided with the -o option. There you will find a collection of assemblies performed by the different programs. The ones named ragoutXXXX are the combined assemblies and it is up to you to choose which assembly is the one you prefer to use. Ragout will only join assemblies that have a similar genome size to the one provided. If none of the assemblies reach this size, Ragout will not run and you will be left with a set of primary assemblies.

