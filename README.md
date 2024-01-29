# skandiver: a divergence-based analysis tool for identifying intercellular mobile genetic elements
## Introduction
skandiver is a program for identifying mobile genetic elements (prophages, plasmids, transposases, etc.) from assembled whole genome sequences using average nucleotide identity (ANI) and evolutionary divergence time.  

## Requirements 
1. skani (Version 0.2.1 or higher)
2. Database of representative genomes (Recommended GTDB Database, see below for installation details)
3. Python 3 (with pandas and bio packages)

## Setting up skani
skandiver uses skani (Developed by Jim Shaw at https://github.com/bluenote-1577/skani), a scalable and robust search tool for computing average nucleotide identity between whole genomes. Ensure that skani is properly installed and in PATH. 
