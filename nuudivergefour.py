#Evolutionary Divergence Annotator V0.04 

import sys
from Bio import Phylo
from Bio import Entrez, SeqIO
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
    #Split the name into words and filter
    parts = [part for part in name.split() if not any(c.isdigit() or c in '-_:;!`~+' for c in part)]
    processed_name = '_'.join(parts[:2]) 
    return processed_name
    
def get_gene_name(genome_id):
    try:
        # Search NCBI database for gene name based on genome ID
        handle = Entrez.esearch(db="nucleotide", term=genome_id)
        #print("handle", handle)
        record = Entrez.read(handle)
        #print("record", record)
        if int(record['Count']) > 0:
            # Fetch the first record (usually the most relevant)
            gene_id = record['IdList'][0]
            summary = Entrez.esummary(db="nucleotide", id=gene_id)
            summary_record = Entrez.read(summary)
            gene_name = summary_record[0]['Title']
            return gene_name
        else:
            return None
    except Exception as e:
        print(f"Error occurred while fetching gene name for {genome_id}: {e}")
        return None

def update_titles(input_file, output_file):
    # Read multifasta file and update titles
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            genome_id = record.id
            print("genome_id", genome_id)
            gene_name = get_gene_name(genome_id)
            if gene_name:
                record.description = gene_name
            SeqIO.write(record, outfile, "fasta")


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
                #In this case, we are unable to determine the divergence time 
                #Possibly due to uncharacterized reference or query species
                #Set divergence time to negative 1 for later post-processing
                divergence_time = -1.0
        else:
            #In this case the reference and query species are the same, thus divergence time is 0
            divergence_time = 0.0

        if isinstance(divergence_time, (float, int)) and divergence_time != 0.0:
            outfile.write(line.strip() + f'\t{divergence_time}\n')

        #Periodically update the progress bar
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
def analyze_results(input_file, output_file, tree, species_file):
    #Read the species names from the text file
    with open(species_file, 'r') as species_file:
        species_names = set(line.strip() for line in species_file)
    df = pd.read_csv(input_file, sep='\t')

    #Store information for each query
    query_info = {}

    #Iterate through dataframe rows
    for index, row in df.iterrows():
        query_name = row['Query_name']
        divergence_time = row['DivergenceTime(MYA)']
        match_name = process_name(row['Ref_name'])
        ani = float(row['ANI'])
        align_frac_ref = float(row['Align_fraction_ref'])
        
        #IMPORTANT: Check if species exist in directory
        if query_name in query_info:
            query_info[query_name]['NumberHits'] += 1
            query_info[query_name]['TotalDivergence'] += divergence_time
            query_info[query_name]['TotalANI'] += ani
            query_info[query_name]['TotalAlignFracRef'] += align_frac_ref
            #Append total species matches 
            if match_name not in query_info[query_name]['RefSpeciesHits']:
                query_info[query_name]['RefSpeciesHits'].add(match_name)
        else:
            query_info[query_name] = {'NumberHits': 1, 'TotalDivergence': divergence_time, 'RefSpeciesHits': {match_name}, 'TotalANI': ani, 'TotalAlignFracRef': align_frac_ref}

    #Recalculate divergence time for unknown query
    for query_name, info in query_info.items():
        if info['TotalDivergence'] < 0.0:
            total_divergence = 0.0
            total_pairs = 0
            
            unknown_species_hits = list(info['RefSpeciesHits'])
            
            #Find the first species in unknown_species_hits that is also in species_names
            selected_species = next((species for species in unknown_species_hits if species in species_names), None)
            
            if selected_species is not None:
                #Calculate pairwise divergence time between selected_species and other species in unknown_species_hits
                for species in unknown_species_hits:
                    if species != selected_species:
                        temporary_divergence_time = get_divergence_time(tree, selected_species, species)
                        
                        if temporary_divergence_time is not None:
                            total_divergence += temporary_divergence_time
                            total_pairs += 1
            
            
            #DEPRECATED PAIRWISE DIVERGENCE
            #Calculate pairwise divergence time for all valid species
            """
            for i in range(len(unknown_species_hits)):
                for j in range(i + 1, len(unknown_species_hits)):
                    species1 = unknown_species_hits[i]
                    species2 = unknown_species_hits[j]
                    if species1 in species_names and species2 in species_names:
                        temporary_divergence_time = get_divergence_time(tree, species1, species2)

                        #print("New divergence time is:", temporary_divergence_time)
                    
                        if temporary_divergence_time is not None:
                            total_divergence += temporary_divergence_time
                            total_pairs += 1
            
            
            if total_pairs > 0:
                new_divergence_time = total_divergence / total_pairs
                info['TotalDivergence'] = new_divergence_time
            """
            
    #Create summary dataframe
    summary_df = pd.DataFrame.from_dict(query_info, orient='index')
    summary_df.reset_index(inplace=True)
    summary_df.columns = ['Query Name', 'NumberHits', 'TotalDivergence', 'RefSpeciesHits', 'TotalANI', 'TotalAlignFracRef']

    #Ensure there exists more than one match for queries
    summary_df = summary_df[summary_df['TotalDivergence'] != -1.0]
    
    #Calculate the average divergence for each query
    #Calculate average ANI/AF
    summary_df['AverageHitDivergence'] = summary_df['TotalDivergence'] / summary_df['NumberHits']
    summary_df['AverageANI'] = summary_df['TotalANI']/summary_df['NumberHits']
    summary_df['AverageAlignFracRef'] = summary_df['TotalAlignFracRef']/summary_df['NumberHits']
    
    #Convert ref hits to string and remove curly brackets and single quotation marks
    summary_df['RefSpeciesHits'] = summary_df['RefSpeciesHits'].apply(lambda x: ', '.join(x))

    split_df = summary_df.apply(extract_query_info, axis=1)
    concat_df = pd.concat([split_df, summary_df], axis=1)
    concat_df = concat_df.drop('Query Name', axis=1)
    
    #Reorder results columns
    output_order = ['GenomeID/AccessionNumber', 'QuerySpecies', 'GenomePosition', 'NumberHits', 'TotalDivergence', 'AverageHitDivergence', 'AverageANI', 'AverageAlignFracRef', 'RefSpeciesHits'] 
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

    analyze_results(output_results_file, output_summary_file, tree, "gtdb_species_addenumall.txt")
    print("Done")
