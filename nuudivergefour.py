#Evolutionary Divergence Annotator V0.04 

import sys
from Bio import Phylo
import re
import pandas as pd

#Function to extract GTDB species names from a column
#DEPRECATED
def extract_species_name(column):
    match = re.search(r'\((\S+ \S+)', column)
    if match:
        return match.group(1).replace(" ", "_")
    else:
        return None

#Function to read in Newick file representing a divergence time tree
def parse_newick_file(newick_file):
    return Phylo.read(newick_file, 'newick')

#Function to navigate the divergence time tree and retrieve the evolutionary distance between two species
def get_divergence_time(tree, species1, species2):
    if not (tree.find_any(species1) and tree.find_any(species2)):
        return None
    return tree.distance(species1, species2) / 2
	
#Function to process a name by removing any numbers/special characters
def process_name(name):
    #Split the name into words and filter out unwanted words
    parts = [part for part in name.split() if not any(c.isdigit() or c in '-_:;!`~+' for c in part)]
    processed_name = '_'.join(parts[:2]) 
    return processed_name

def add_divergence_time_column(tree, lines, outfile, species_file):
    #Read the species names from the text file
    with open(species_file, 'r') as species_file:
        species_names = set(line.strip() for line in species_file)
    
    #Write the header line with an additional column
    outfile.write(lines[0].strip() + '\tDivergenceTime(MYA)\n')
    total_lines = len(lines)-1

    for idx, line in enumerate(lines[1:], start=1):
        cols = line.split('\t')
        reference_species = process_name(cols[5])
        query_species = process_name(cols[6])

        if reference_species != query_species:
            if reference_species in species_names and query_species in species_names:
                divergence_time = get_divergence_time(tree, reference_species, query_species)
            else:
                divergence_time = 0.0
        else:
            divergence_time = 0.0

        if isinstance(divergence_time, (float, int)) and divergence_time != 0.0:
            outfile.write(line.strip() + f'\t{divergence_time}\n')

        # Update the custom progress bar
        progress = int(idx / total_lines * 100)
        bar_length = 50
        bar = '[' + '>' * int(bar_length * idx / total_lines) + ' ' * (bar_length - int(bar_length * idx / total_lines)) + ']'
        sys.stdout.write(f'\rProgress: {progress}% {bar}')
        sys.stdout.flush()

def extract_query_info(row):
    # Example function to extract information from the Query Name
    query_parts = row['Query Name'].split()

    # Assuming that GenomeID/AccessionNumber is the first part before the first space
    genome_id = query_parts[0]

    # Assuming that QuerySpecies is the parts between the second and last words
    #query_species = ' '.join(query_parts[1:-6])
    query_species = process_name(row['Query Name'])

    # Assuming that Genome Position is the portion after the last colon
    # and excluding any leading or trailing characters such as parentheses
    genome_position = query_parts[-1].split(':')[-1].strip('()')

    return pd.Series({'GenomeID/AccessionNumber': genome_id,
                      'QuerySpecies': query_species,
                      'GenomePosition': genome_position})

#Function to output a set of potential MEs (fragments with significant total divergence time)
def analyze_results(input_file, output_file):
    df = pd.read_csv(input_file, sep='\t')

    #Store information for each query
    query_info = {}

    #Iterate through dataframe rows
    for index, row in df.iterrows():
        query_name = row['Query_name']
        divergence_time = row['DivergenceTime(MYA)']
        match_name = process_name(row['Ref_name'])
        #IMPORTANT: Check if species exist in directory
        if query_name in query_info:
            query_info[query_name]['NumberHits'] += 1
            query_info[query_name]['TotalDivergence'] += divergence_time
            if match_name not in query_info[query_name]['RefSpeciesHits']:
                query_info[query_name]['RefSpeciesHits'].add(match_name)
        else:
            query_info[query_name] = {'NumberHits': 1, 'TotalDivergence': divergence_time, 'RefSpeciesHits': {match_name}}

    #Create summary dataframe
    summary_df = pd.DataFrame.from_dict(query_info, orient='index')
    summary_df.reset_index(inplace=True)
    summary_df.columns = ['Query Name', 'NumberHits', 'TotalDivergence', 'RefSpeciesHits']

    #Calculate the average divergence for each query
    summary_df['AverageHitDivergence'] = summary_df['TotalDivergence'] / summary_df['NumberHits']
    
    #Convert ref hits to string and remove curly brackets and single quotation marks
    summary_df['RefSpeciesHits'] = summary_df['RefSpeciesHits'].apply(lambda x: ', '.join(x))

    split_df = summary_df.apply(extract_query_info, axis=1)
    concat_df = pd.concat([split_df, summary_df], axis=1)
    concat_df = concat_df.drop('Query Name', axis=1)
    
    #Reorder results columns
    output_order = ['GenomeID/AccessionNumber', 'QuerySpecies', 'GenomePosition', 'NumberHits', 'TotalDivergence', 'AverageHitDivergence', 'RefSpeciesHits'] 
    concat_df = concat_df.reindex(columns=output_order)
    
    #Write the summary DataFrame to a new results file (tabs as separators)
    with open(output_file, 'w') as outfile:
        outfile.write("skandiver Potential Mobile Elements Results\n\n")
        concat_df.to_csv(outfile, sep='\t', index=False, lineterminator='\n')


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python nuudivergefour.py input_file output_file")
        sys.exit(1)

    input_results_file = sys.argv[1]
    output_results_file = "skanichunkingext.txt"
    newick_file_path = "gtdb_species_addenum.nwk"
    output_summary_file = sys.argv[2]

    tree = parse_newick_file(newick_file_path)

    #Read the existing results file
    with open(input_results_file, 'r') as infile:
        lines = infile.readlines()

    #Open the output file for writing
    with open(output_results_file, 'w') as outfile:
        add_divergence_time_column(tree, lines, outfile, "gtdb_species_addenumall.txt")

    analyze_results(output_results_file, output_summary_file)
    print("Done")