# skandiver: a divergence-based analysis tool for identifying intercellular mobile genetic elements
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
*skani is installed and in PATH (i.e. typing ```skani -h``` works). Visit https://github.com/bluenote-1577/skani for more information on setting up skani.
*~120 GB total free disk space for the uncompressed database and indexing.

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


## Installation and Testing 
First download the skandiver repository: 
```sh
git clone https://github.com/YoukaiFromAccounting/skandiver
cd skandiver
bash SETUP.sh
```

The provided setup script will test your environment for dependencies and download an example data set. You can also install all needed dependencies using the following: 
```sh
sudo apt-get install python3-pip
pip3 install bio pandas
```
