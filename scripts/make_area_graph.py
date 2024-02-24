# This code will take in the output .txt file of skandiver
# and output a graph that shows all regions of detected 
# mobile elements in a genome
# must have matplotlib installed

import sys
import pandas as pd
import csv
import os
import subprocess
try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
except ModuleNotFoundError:
    print("Error: Matplotlib is not installed.")
    print("Please install Matplotlib by running: pip install matplotlib")
    exit()
import math

def graph_genome(name_tuple, coordinates, row, length):
    ax = fig.add_subplot(length, 1, row) #**
    #print("running graph_genome")
    ax.set_title(f'Areas of detected mobile elements in {name_tuple[0]} {name_tuple[1]}')
    skani_data = coordinates

    # Create a figure and axis
    #fig, ax = plt.subplots(figsize=(9, 0.5))
    
    
    max_coor = 0
    for start1, end1 in skani_data:
        #print("here", start1, end1)
        rect1 = patches.Rectangle((start1, 0), end1 - start1, 1, linewidth=1, edgecolor= 'red', facecolor='red', alpha=0.5)
        ax.add_patch(rect1)
        if end1 > max_coor:
            max_coor = end1 ##change this to the size of the actal genome eventually

    # Hide y-axis
    ax.set_yticks([])
    ax.set_yticklabels([])

    # Add labels
    #print("max_coor", max_coor)
    step = max_coor//10
    boundary = max_coor+step
    ax.set_xticks(range(0,  boundary, step))
    ax.set_xlim(0, boundary)
    #print("xticks set")
    ax.set_xlabel('Genomic Position')
    #print("xlabel set")
    
    #axs[row].set_ylim()
    


    

def get_coor(file_path):
    #print("running get_coor")

    # dictionary { (genomeID, query_species):[list of coordinates] }
    genome_dict = {}
    # Open the file for reading
    with open(file_path, 'r') as file:
        # Initialize an empty list to store the extracted genome positions
        
        # Iterate over each line in the file
        for line in file:
            #print("LINE", line)
            if len(line) < 3: continue
            # Split the line into columns based on whitespace
            columns = line.split()
            
            # Extract the GenomePosition column (assuming it's the third column)
            genome_position = columns[2]
            genomeID = columns[0]
            query_species = columns[1]
            #print("genome_postion", genome_position)
            if genome_position == 'Element' or genome_position == 'GenomePosition': continue
            if genome_position == 'Mobile': continue
            #print("Here")
            id_tuple = (genomeID, query_species)
            if id_tuple not in genome_dict:
                genome_dict[id_tuple] = []
            
            # Convert the genome position string into a tuple
            start, end = map(int, genome_position.split('-'))
            #print("start, end", start, end)
            genome_dict[id_tuple].append((start, end))
    
    # Output the list of genome positions
    return genome_dict


if __name__ == "__main__":
    # Check if the correct number of command-line arguments is provided
    if len(sys.argv) != 2:
        print("Usage: python3 make_area_graph.py <path_to_skandiver.txt_ouput>")
        sys.exit(1)

    # Get command-line arguments
    file_path = sys.argv[1]
    file_name = os.path.basename(file_path)

    genome_positions = get_coor(file_path)
    #print("extract genome_positions done")
    #print("genome_positions: ", genome_positions)

    num_genes = len(genome_positions)
    num_rows = 6
    if num_genes <= num_rows:
        fig = plt.figure(figsize=(9, 1.5*num_genes))  

        row = 1
        for gene in genome_positions:
            graph_genome(gene, genome_positions[gene], row, len(genome_positions))
            row += 1

        plt.tight_layout()
        #plt.subplots_adjust(hspace=2)
        plt.savefig(f'{file_name}_plot.png')
        plt.show()
    else:
        print("We can only accommodate 6 graphs at a time on one plot. Each plot will be saved in separate files.")
        key_vals = list(genome_positions.items())
        num_figs = math.ceil(num_genes / num_rows)
        gene_num = 0
        for i in range(num_figs):
            fig = plt.figure(figsize=(9, 1.5*num_rows))  

            row = 0
            for gene_index in range(gene_num, min(gene_num+num_rows, num_genes)):
                gene_name = key_vals[gene_index][0]
                gene_vals = key_vals[gene_index][1]
                graph_genome(gene_name, gene_vals, row+1, num_rows)
                row = (row+1)%num_rows

            plt.tight_layout()
            #plt.subplots_adjust(hspace=2)
            plt.savefig(f'{file_name}_plot.png')
            plt.show()