# Phylogenomics of PERU and FLS2 LRR-RK immune receptors in Solanaceae plant family
# Part one - Database preparation

After downloading the data from [proteomes](/proteomes) unzip it using the following command (mac OS) or other decompressing programs:

```bash 
cd /path/to/proteomes/ 

gunzip *.gz
```

Then in R compile all into a single protein fasta file, which will be used to extract LRR-RK sequences from:

```R
# R

# load the requierd packages
library(Biostrings)
library(tidyverse)
library(readr)

# set the path to the directory containing the proteome FASTA files
fasta_dir <- "/path/to/the/proteoms/"

# get a list of all the file names in the directory
fasta_files <- list.files(fasta_dir, pattern = ".fasta", recursive = TRUE)

# use lapply() to read in each file and store the sequences as a DNAStringSet object
fasta_seqs <- lapply(fasta_files, function(x) {
  # Read in the fasta file
  seqs <- readAAStringSet(file.path(fasta_dir, x))
  
  # Add the file name as metadata to each sequence
  mcols(seqs)$file_name <- x
  
  # Return the AAStringSet object
  seqs
})

# combine all the sequences into a single AAStringSet object
SolDB <- do.call(c, fasta_seqs)
SolDB@ranges@NAMES <- sub(" .*", "", SolDB@ranges@NAMES) # protein name clean-up

# metadata dataframe
SolDB_meta <- data.frame(ID = SolDB@ranges@NAMES,
                         genome = SolDB@elementMetadata@listData)
SolDB_meta$file_name <- gsub("*.fasta", "", SolDB_meta$file_name) # file name clean-up

# write the DB
writeXStringSet(SolDB, "/path/to/your/directory/SolDB.fasta")
```

Import the genome data on [table S2](/tables/table_S2.csv) and compile the metadata:

```R
# R

# metadata compilation

# import the genome metadata information file (table S2)
Genome_data <- read_csv("/path/to/tables/table_S2.csv") # import genome metadata

genome_meta <- Genome_data[,c(1,  # organism
                              2,  # scientific name
                              3)]  # file_name

# merge the database metadata and genome metadata based on file names to have everything all together
SolDB_metadata <- SolDB_meta %>% left_join(genome_meta, by = "file_name")

write_csv(SolDB_metadata, "/path/to/your/directory/SolDB_metadata.csv") # export if you like
```

[Main](README.md) | [Next: 02. Search for PK and LRR containing proteins](02_Search_for_PK_and_LRR_containing_proteins.md)
