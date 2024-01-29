# skandiver: a divergence-based analysis tool for identifying intercellular mobile genetic elements
# Current version: v0.1.1 - 2024-01-25
## Introduction
skandiver is a program for identifying mobile genetic elements (prophages, plasmids, transposases, etc.) from assembled whole genome sequences using average nucleotide identity (ANI) and evolutionary divergence time. 

## Requirements 
1. skani (Version 0.2.1 or higher)
2. Database of representative genomes (Recommended GTDB Database, see below for installation details)
3. Python 3 (with pandas and bio packages)

## Setting up skani
skandiver uses skani (Developed by Jim Shaw at https://github.com/bluenote-1577/skani), a scalable and robust search tool for computing average nucleotide identity between whole genomes. Ensure that skani is properly installed and in PATH. 

## Setting up database of representative genomes
skani search requires a database of representative genomes to query against. The current recommended database is the Genome Taxonomy Database (GTDB), which contains >85,000 representative genomes. 
To setup this database, first ensure that the following requirements are met: 
* skani is installed and in PATH (i.e. typing ```skani -h``` works). Visit https://github.com/bluenote-1577/skani for more information on setting up skani.
* ~120 GB free disk space is available for the uncompressed database and indexing.

First, download the compressed GTDB database and unzip it:
```
wget https://data.gtdb.ecogenomic.org/releases/release214/214.1/genomic_files_reps/gtdb_genomes_reps_r214.tar.gz
tar -xf gtdb_genomes_reps_r214.tar.gz
```

The gtdb database is formatted in a special way. In order to process the reference genome files inside the gtdb folder, we have to do a bit of work. We can run the following to collect all genomes locations into a file called ```gtdb_file_names.txt```.
```
find gtdb_genomes_reps_r214/ | grep .fna > gtdb_file_names.txt
```

Next, we can construct the indexed database to query against using:
```
skani sketch -l gtdb_file_names.txt -o gtdb_skani_database_ani -t 20
```


## Installation and Quick Start
First download the skandiver repository: 
```sh
git clone https://github.com/YoukaiFromAccounting/skandiver
cd skandiver
chmod +x skandiver.sh
bash SETUP.sh
```

The provided setup script will test your environment for dependencies and download an example data set. You can also install all needed dependencies using the following: 
```sh
sudo apt-get install python3-pip
pip3 install bio pandas
```

You can test skandiver against a sample whole genome assembly of _Acinetobacter baumannii_ by executing the following command: 
```
./skandiver.sh test_files results 10000 [PATH_TO_REPRESENTATIVE_GENOME_DB]
```
For example, if you followed the above instructions for setting up the GTDB database of representative genomes in the skandiver directory, you can run:
```
./skandiver.sh test_files results 10000 gtdb_skani_database_ani
```
This should output four files; results.txt, resultsskani.txt, resultsskanifiltered.txt, and resultssearch.fna. results.txt contains the summary of potential mobile genetic elements found by skandiver, while resultsskani.txt and resultsskanifiltered.txt contain the skani search results for the query whole genome assembly (with resultsskanifiltered only displaying genome matches with greater than 95% average nucleotide identity and 90% align fraction). resultssearch.fna contains the entire fragmented genome assembly used for the skani search. 

## Contact
  Brian Zhang, xiaoleiz@andrew.cmu.edu (Contributing author)  
  Grace Oualline, gouallin@andrew.cmu.edu (Contributing author)

## Acknowledgements
We would like to express our gratitude to the following individuals and organizations for their major contributions and support in the development of skandiver: 
* Jim Shaw (https://github.com/bluenote-1577) for the creation and continuous support of skani, a fundamental tool utilized by skandiver for ANI computations, as well as providing valuable guidance regarding the overall quality and usability of skandiver.
* Yun William Yu (https://github.com/yunwilliamyu) for providing algorithmic support and troubleshooting expertise, greatly improving skandiver's efficiency.

This implementation of skandiver was based on the ideas and software from the following paper:   
    Shaw, J., & Yu, Y. W. (2023). Fast and robust metagenomic sequence comparison through sparse chaining with skani. Nature methods, 20(11), 1661–1665. https://doi.org/10.1038/s41592-023-02018-3
