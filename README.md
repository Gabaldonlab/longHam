# longHam

LongHam is an assembly pipeline designed to assemble small fungal genomes based on a combination of long and short reads. The pipeline proceeds as follows:

1.- Read cleaning and filtering

2.- Construction of primary assemblies

3.- Fusion of primary assemblies

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
