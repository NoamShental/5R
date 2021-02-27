5R software package![GitHub Logo](logo.png)
========================
**5R package** – A software package accompanying a novel method of bacterial 16S rRNA profiling based on the 5R protocol. Five short regions along the 16S rRNA are amplified in multiplex and subsequently sequenced by an Illumina machine. Resulting reads from these regions are combined into a coherent solution using our Short MUltiple Regions Framework (SMURF) algorithm (Fuks et al. 2018).

Two options are provided: A Matlab code and a standalone executable version for Linux.
The input to either methods is the same, and hence described below for the Matlab version followed by a short description of the standalone syntax.

Installation
------------
Download the 5R package.


Usage
-------
Open Matlab, change to the mFiles directory and run the following line:

main_5R(illumina_files_dir, db_dir, db_name, results_filename, kmer_len, generate_group_results)

where the inputs are:

**illumina_files_dir** – a (**full**) path to directory of fastq paired-end sequencing files of samples profiled by the 5R primers. The package supports two cases: (i) the fastq files  corresponding to each sample are in a different sub-directory or (ii) fastq of different samples are mixed together in one folder. In the latter case the package will automatically split the files to folders per sample based on the the file name.

**db_dir** - the (**full**) path to the directory of data base kmers

**db_name** - the name of the data base (the default data base of 5R is name is "GreenGenes_201305")

**results_filename** - a **full** path to the results file name. Remark: Note that the directory in which results would be written must be created in advance.

**kmer_len** - the length of the k-mer to be applied for reconstruction. The maximal supported length is 160nt. To support higher resolution the kmer length should be set as high as possible, yet depending on read quality. Since read qualify often deteriorate with the reads' length, we suggest the length over which the typical read quality falls below 30.
If not specified, the a default of 100nt is set.

**generate_group_results** - (this input parameter is optional) - this is a 0/1 indicator. When true, a folder named "groups" will be created inside the "resDir" folder for each sample and the full 16S sequences of all the reconstructed groups (i.e. amplicons) will be written to fasta files. A separate file will be created for each group containing all the sequences of the group.



Output
-------
Two files are created:

**results_filename.txt** - A table of species' relative abundances (rows) in each sample (columns).

**ReadCountStats_results_filename.txt** - General statistics of the number reads that passed the quality filter, the number of reads matched to each regions etc.


Example
-------

The directory "example_fastq" contains 5R sequencing files of two breast cancer FFPE samples.
To perform SMURF reconstruction of these samples use the following command

main_5R('./example_fastq','./GG_5R', './example_results/5R_SMURF_example.txt',126);  


Standalone vesion
----------------
The 5R package can also be applied on a linux machine without a Matlab license.

Download the following MCR file from the Mathworks site:

https://ssd.mathworks.com/supportfiles/downloads/R2019b/Release/4/deployment_files/installer/complete/glnxa64/MATLAB_Runtime_R2019b_Update_4_glnxa64.zip

Following download double click on the installer and follow the instructions. This package is required for running standalone code compiled using Matlab.

For instructions regarding the MCR installation see:

https://www.mathworks.com/help/compiler/install-the-matlab-runtime.html




To run the 5R code use the following syntax

./5R_linux/run_main_5R.sh PATH_TO_RUNTIME illumina_files_dir db_dir  db_name  results_filename <kmer_len>  <1>

where PATH_TO_RUNTIME is the directory in which the MCR was installed (e.g. /opt/apps/matlab/Runtime/v97). Replace PATH_TO_RUNTIME by your directory.

For example, the following syntax would work on the two above mentioned samples:

./5R_linux/run_main_5R.sh PATH_TO_RUNTIME  ./example_fastq  ./GG_5R GreenGenes_201305 ./example_results/5R_SMURF_example.txt 126


Contact us
----------------
For questions please email:

Garold Fuks: garoldf@gmail.com

Noam Shental: shental@openu.ac.il

Ravid Straussman: ravidst@weizmann.ac.il



## References

Garold Fuks, Michael Elgart, Amnon Amir, Amit Zeisel, Peter J. Turnbaugh, Yoav Soen and Noam Shental. Combining 16S rRNA gene variable regions enables high-resolution microbial community profiling. *Microbiome*, 26;6(1):17, Jan 2018. doi: 10.1186/s40168-017-0396-x.

