# Standardize_cpDNA
### Standardization of the chloroplast DNA or Plastome assemblies

Standardizing the cpDNA assemblies makes it easier to run many downstream analyses, especially plastome comparisons. This simple script is written to standardize any cpDNA assembly to start from an LSC region followed by IRa, SSC, and IRb regions. 

### *Pre-requisites* 

Install the below software and add them to the environment

> samtools (tested with v1.13)

> blast+ (tested with v2.12.0+)

### Install & run

    git clone https://github.com/SaiReddy-A/Standardize_cpDNA.git
    
    ./standardize_cpDNA.sh -i input.fa
    
The output is a standardized fasta file written to `input.standardized.fa`

### Example 

    ./standardize_cpDNA.sh -i test.fa
    


<br />

<br />

<br />



*Any suggessions and comments are greatly appreciated!!*

### *Thank you!*
