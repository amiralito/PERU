After compiling the protein database, extract the LRR-RK using the following steps:

1. Search for proteins containing PK domain using *hmmsearch* and PF00069.hmm
```bash
# first search for Protein Kinase domain by searching the dataset for Pfam PF00069 hmm profile using hmmsearch and following options
hmmsearch -o PF00069.out --tblout PF00069_tbl.out --domtblout PF00069_domtbl.out --max PF00069.hmm SolDB.fasta
```

- Import and process the *hmmsearch* output:
```R
library(Biostrings)
library(tidyverse)
library(readr)

# import the PK hmmsearch output
PK_tbl <- read_table("/path/to/PF00069_tbl.out",
                     col_names = FALSE, col_types = cols(X19 = col_skip()), 
                     skip = 3) %>% head(-10)


PK_scores <- PK_tbl[,c(1,5,6,7)] %>% setNames(c("ID","E-Value","score","bias")) # dataframe with only scores and e-values
PK_scores_filtered <- PK_scores[PK_scores$`E-Value` < 0.01,] # remove anything below the e-value threshold of 0.01

PK_seq <- SolDB_geneious[PK_scores_filtered$ID] # PK containing sequences


# keep onky the sequences with one kinase domain
# import the domain table output from hmmsearch
PK_domtbl <- read_table("/path/to/PF00069_domtbl.out",
                        col_names = FALSE, col_types = cols(X19 = col_skip()), 
                        skip = 3) %>% head(-10)



PK_domtbl_filtered <- filter(PK_domtbl, PK_domtbl$X11 == 1) # keep only the ones with one PK domain
PK_domtbl_filtered <- filter(PK_domtbl_filtered, PK_domtbl_filtered$X7 < 0.01) # remove anything below the e-value threshold of 0.01

# refine the PK table and export it for next step
PK_filtered <- PK_scores_filtered[PK_scores_filtered$ID %in% PK_domtbl_filtered$X1,]
PK_filtered_seq <- PK_seq[PK_filtered$ID]

write_csv(PK_filtered, "/path/to/PK_filtered.csv")
writeXStringSet(PK_filtered_seq, "/path/to/PK_filtered.fasta")

# extract the PK domains only. This will be used for phylogenetic analysis
PK_domain_seq <- subseq(PK_filtered_seq, start = PK_domtbl_filtered$X20, end = PK_domtbl_filtered$X21)

```

2. Search for presence of LRR domain in the PK containing sequences using the PF00560.hmm, PF08263.hmm, PF13516.hmm, and PF13855.hmm
```bash
# LRR_1 (Pfam PF00560)
hmmsearch -o PF00560.out --tblout PF00560_tbl.out --domtblout PF00560_domtbl.out --max PF00560.hmm PK_filtered.fasta

# LRRNT_2 (Pfam PF08263)
hmmsearch -o PF08263.out --tblout PF08263_tbl.out --domtblout PF08263_domtbl.out --max PF08263.hmm PK_filtered.fasta

# LRR_6 (Pfam PF13516)
hmmsearch -o PF13516.out --tblout PF13516_tbl.out --domtblout PF13516_domtbl.out --max PF13516.hmm PK_filtered.fasta

# LRR_8 (Pfam PF13855)
hmmsearch -o PF13855.out --tblout PF13855_tbl.out --domtblout PF13855_domtbl.out --max PF13855.hmm PK_filtered.fasta
```

- Import and process the *hmmsearch* output
```R
library(Biostrings)
library(tidyverse)
library(readr)

# LRR1
LRR1_tbl <- read_table("/path/to/PF00560_tbl.out",
                        col_names = FALSE, col_types = cols(X19 = col_skip()), 
                        skip = 3) %>% head(-10)
LRR1_tbl <- filter(LRR1_tbl, LRR1_tbl$X5 < 0.01) # remove anything below the e-value threshold of 0.01
LRR1 <- LRR1_tbl$X1 %>% as.data.frame() # keep only the ids

#LRR6
LRR6_tbl <- read_table("/path/to/PF13516_tbl.out",
                       col_names = FALSE, col_types = cols(X19 = col_skip()), 
                       skip = 3) %>% head(-10)
LRR6_tbl <- filter(LRR6_tbl, LRR6_tbl$X5 < 0.01) # remove anything below the e-value threshold of 0.01
LRR6 <- LRR6_tbl$X1 %>% as.data.frame() # keep only the ids

# LRR8
LRR8_tbl <- read_table("/path/to/PF13855_tbl.out",
                       col_names = FALSE, col_types = cols(X19 = col_skip()), 
                       skip = 3) %>% head(-10)
LRR8_tbl <- filter(LRR8_tbl, LRR8_tbl$X5 < 0.01) # remove anything below the e-value threshold of 0.01
LRR8 <- LRR8_tbl$X1 %>% as.data.frame() # keep only the ids

#LRRNT2
LRRNT2_tbl <- read_table("/path/to/PF08263_tbl.out",
                      col_names = FALSE, col_types = cols(X19 = col_skip()), 
                      skip = 3) %>% head(-10)
LRRNT2_tbl <- filter(LRRNT2_tbl, LRRNT2_tbl$X5 < 0.01) # remove anything below the e-value threshold of 0.01
LRRNT2 <- LRRNT2_tbl$X1 %>% as.data.frame() # keep only the ids

# merge all LRR containing proteins
LRR <- Reduce(function(x,y) merge(x,y, by = ".", all = TRUE), list(LRR1,LRR6,LRR8,LRRNT2)) # merge all LRR containing proteins

names(LRR) <- "ID"

# subset the PK domain of the remaining proteins for alignment and phylogenetic analysis
LRR_seq <- PK_filtered_seq[LRR$ID]
LRR_PK_seq <- PK_domain_seq[LRR$ID]

writeXStringSet(LRR_seq, "/path/to/LRR.fasta") # full-length
writeXStringSet(LRR_PK_seq, "/path/to/LRR_PKD.fasta") # PK domain only
```