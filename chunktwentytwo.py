#Short chunking script V0.22 (added directory functionality) 

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import os
import glob
from tqdm import tqdm

#Take in user input for algorithm parameters
if len(sys.argv) != 4:
    print("Incorrect script call, ensure call is formatted as: python3 chunketc.py <input_directory> <output directory> <chunksize>")
    sys.exit(1)
  
#Read in name of the file and the desired size of chunks in nucleotides
input_directory = sys.argv[1]
output_file = sys.argv[2]
chunk_size = int(sys.argv[3])

#Use SeqIO.parse instead of read for multifasta file IMPORTANT
#records = list(SeqIO.parse(input_file, "fasta")) <<< moved to multi_batch function 

def create_batch(records, chunk_size):
    record_it = iter(records)

    record = next(record_it)
    current_base = 0

    batch = []
    batch_size = 0

    while record:
        #Loop over records untill the batch is full. (or no new records)
        while batch_size != chunk_size and record:

            end = current_base + chunk_size - batch_size
            seq = record[current_base:end]

            end_of_slice = current_base + len(seq) - 1
            #Use record.id for pure name, record.description for contig details
            fasta_header = record.description + ":{}-{}".format(current_base, end_of_slice)
            
            seq.id = seq.name = fasta_header
            seq.description = ''
            batch.append(seq)

            current_base += len(seq)
            batch_size += len(seq)

            #Current record is exhausted, get a new one.
            if current_base >= len(record):
                record = next(record_it, None)
                current_base = 0

        #Extracted batch with the correct size
        yield batch
        batch = []
        batch_size = 0
        
        
def group_batch(input_directory, output_file, chunk_size):
    #Ensure output directory exists
    output_dir = os.path.dirname(output_file)
    
    if output_dir:
        os.makedirs(output_dir, exist_ok = True)
    
    #Initialize a list to store chunks
    all_chunks = []
    
    #Obtain total files to process
    total_files = len(glob.glob(os.path.join(input_directory, '*.fna'))) + len(glob.glob(os.path.join(input_directory, '*.fasta')))
    with tqdm(total=total_files, desc="Processing Files", unit="file") as pbar:
        #Iterate over .fna files individually within input directory 
        for ind_file in glob.glob(os.path.join(input_directory, '*.fna')):
            records = list(SeqIO.parse(ind_file, "fasta"))
            
            #Iterate over individual batches and run chunking 
            for i, batch in enumerate(create_batch(records, chunk_size)):
                all_chunks.extend(batch)
                pbar.update(1)
                #filename = os.path.join(output_directory, f"{os.path.splitext(os.path.basename(ind_file))[0]}_chunk{i}.fasta")
                #SeqIO.write(batch, filename, "fasta") 
                #print("Successfully chunked",batch, filename)
               
        #Iterate over .fasta files individually within input directory 
        for ind_file in glob.glob(os.path.join(input_directory, '*.fasta')):
            records = list(SeqIO.parse(ind_file, "fasta"))
            
            #Iterate over individual batches and run chunking 
            for i, batch in enumerate(create_batch(records, chunk_size)):
                all_chunks.extend(batch)
                pbar.update(1)
                #filename = os.path.join(output_directory, f"{os.path.splitext(os.path.basename(ind_file))[0]}_chunk{i}.fasta")
                #SeqIO.write(batch, filename, "fasta") 
                #print("Successfully chunked",batch, filename)
    
    #Write chunks to singular multifasta file
    SeqIO.write(all_chunks, output_file, "fasta")
    print("Successfully created multifasta file:", output_file)


group_batch(input_directory, output_file, chunk_size)

#DEPRECATED
#Ensure output directory exists
#os.makedirs(output_directory, exist_ok = True)

#Range over batches and export each chunk as a fasta file
#for i, batch in enumerate(create_batch(records, chunk_size)):
    #Add input_file to prevent file overwrite from chunk loading
    #filename = os.path.join(output_directory,(input_file+"chunk{}.fasta").format(i))
    #SeqIO.write(batch, filename, "fasta")