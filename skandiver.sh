#!/bin/bash

#Function ensuring files/directories are present
check_directory() {
    if [ ! -d "$1" ]; then
        echo "Error: $2 not found or incorrectly formatted; please double-check if file initialized correctly."
        echo "Exiting skandiver."
        exit 1
    fi
}

#Function to extract fasta/fna files from compressed files
gunzip_files() {
    echo "Step 1: Gunzipping files..."
    for file in "$1"/*.gz; do
        gunzip "$file"
    done
}

#Function to fragment genomes according to a set size
fragment_genomes() {
    echo "Step 2: Fragmenting genomes..."
    echo "python3 chunktwentytwo.py $1 $2 $3"
    python3 chunktwentytwo.py "$1" "$2" "$3"
}

#Function to search fragmented genomes against representative genome database
skani_search() {
    echo "Step 3: skani searching against representative genomes..."
    echo "skani search --qi $1 -d $2 -o $3 -t 10"
    skani search --qi "$1" -d "$2" -o "$3" -t 10
}

#Extract search results with a high match
filter_skani_results() {
    echo "Step 4: Filtering skani search result..."
    echo "awk 'NR==1 || (\$3 > 95 && \$5 > 90)' $1 > $2"
    awk 'NR==1 || ($3 > 95 && $5 > 90)' "$1" > "$2"

    #Check if filtered file can be annotated
    if [ $(wc -l < "$2") -eq 1 ]; then
        echo "No results to annotate in $2. Exiting..."
        exit 1
    fi
}

#Annotate search results with relevant evolutionary divergence information
annotate_divergence() {
    echo "Step 5: Annotating evolutionary divergence information..."
    echo "python3 nuudivergefour.py $1 $2"
    python3 nuudivergefour.py "$1" "$2.txt"
}

main() {
    if [ "$#" -ne 4 ]; then
        echo "Usage: $0 INPUT_DIRECTORY OUTPUT_NAME CHUNKSIZE REFERENCE_DATABASE"
        exit 1
    fi

    INPUT_DIRECTORY=$1
    OUTPUT_NAME=$2
    CHUNKSIZE=$3
    REFERENCE_DATABASE=$4

    OUTPUT_NAME_SEARCH="${OUTPUT_NAME}search.fna"
    OUTPUT_NAME_SKANI="${OUTPUT_NAME}skani.txt"
    OUTPUT_NAME_SKANI_FILTERED="${OUTPUT_NAME}skanifiltered.txt"

    check_directory "$INPUT_DIRECTORY" "Input directory"
    check_directory "$REFERENCE_DATABASE" "Reference database"
	
	echo ========================================
    gunzip_files "$INPUT_DIRECTORY"
	
	echo ========================================
    fragment_genomes "$INPUT_DIRECTORY" "$OUTPUT_NAME_SEARCH" "$CHUNKSIZE"
	
	echo ========================================
    skani_search "$OUTPUT_NAME_SEARCH" "$REFERENCE_DATABASE" "$OUTPUT_NAME_SKANI"
	
	echo ========================================
    filter_skani_results "$OUTPUT_NAME_SKANI" "$OUTPUT_NAME_SKANI_FILTERED"
	
	echo ========================================
    annotate_divergence "$OUTPUT_NAME_SKANI_FILTERED" "$OUTPUT_NAME"
	
	echo ========================================
    if [ -s "${OUTPUT_NAME}.txt" ]; then
        echo "skandiver completed successfully. Potential mobile element regions written to ${OUTPUT_NAME}.txt"
    else
        echo "Error: skandiver was unable to find any mobile regions of interest in $INPUT_DIRECTORY"
    fi
}

main "$@"
