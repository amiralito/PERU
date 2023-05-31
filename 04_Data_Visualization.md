# Phylogenomics of PERU and FLS2 LRR-RK immune receptors in Solanaceae plant family
# Part four - Data Visualization

Once the PERU and FLS2 phylogenetic clades were extracted, heatmaps were generated in R using the following script. The phylogenetic trees were then visualized using iTOL and manually assembled with the heatmaps.

- raw data preparation:
```R
library(Biostrings)
library(tidyverse)
library(readr)

# import the tip labels from the phylogenetic trees made from previous steps
# import the tip labels extracted from the tree in the same order. this was done using TreeViewer software

DM_hm_nr <- read_table("/path/to/trees/PERU.txt", col_names = FALSE) %>% setNames("ID") # PERU clade
FLS2 <- read_table("/path/to/trees/FLS2.txt", col_names = FALSE) %>% setNames("ID") # FLS2 clade


# import the amplified seqeuences metadata (table S6)
hm_meta <- read_csv("/path/to/tables/table S6.csv")
hm_meta$Species <- gsub("S.","Solanum", hm_meta$Species) # slight modifications to make it compatible with the rest of the data


## extract metadata and tidy it up

# PERU and homologs responsive/non-responsive
DM_hm_nr_meta <- SolDB_metadata[SolDB_metadata$ID %in% DM_hm_nr$ID,]
DM_hm_nr_meta <- DM_hm_nr %>% left_join(DM_hm_nr_meta, by = "ID")

DM_hm_nr_meta <- merge(DM_hm_nr_meta, hm_meta, by = "ID", all.x = TRUE) # add the species metadata of the responsive/non-responsive homologs
DM_hm_nr_meta$SciName <- coalesce(DM_hm_nr_meta$SciName.x, DM_hm_nr_meta$SciName.y) # clean up
DM_hm_nr_meta <- DM_hm_nr_meta[,-c(7,8)]

DM_hm_nr_meta <- DM_hm_nr %>% left_join(DM_hm_nr_meta, by = "ID") # reorder the based on the phylogenetic tree


# FLS2
FLS2_meta <- SolDB_metadata[SolDB_metadata$ID %in% FLS2$ID,]
FLS2_meta <- FLS2 %>% left_join(FLS2_meta, by = "ID") # reorder the based on the phylogenetic tree



# import the species list

species_list <- read_csv("/path/to/trees/species_list.txt", col_names = FALSE) %>% 
  setNames("species")

# import the organism list (genomes)

org_list <- read_csv("/path/to/trees/org_list.txt", col_names = FALSE) %>% 
  setNames("Organism")
```

- producing presence/absence matrices

##### Per Species:
This code creates a binary data frame that indicates whether each sequence in the `DM_hm_nr_meta` data frame belongs to a particular species. The [`species_list`](trees/species_list.txt) data frame contains a list of species names. The code iterates over the rows of `species_list` and creates a logical vector that indicates whether the current sequence name matches the species name. The code then creates a column in the binary data frame for the current species name and populates the column with the values from the logical vector. The code then sets the row names of the binary data frame to the sequence names, removes the sequence name column, replaces all NA values with 0, and transposes the data frame. The resulting data frame can be used to identify the species of each sequence in the `DM_hm_nr_meta` data frame. This will be used to generate the presence/absence heatmap. The same code is also used to generate a similar binary data frame for FLS2 clade.

```R
## PERU with responsive and non-responsive homologs

# create a data frame with the sequence names as the only column.
DM_hm_nr_binary_df <- data.frame(seqname = DM_hm_nr_meta$seqname)

# iterate over the rows of the species list.
for (i in 1:nrow(species_list)) {
  species <- species_list$species[i] # get the species name for the current row.
  matches <- ifelse(DM_hm_nr_meta$SciName == species, 1, 0) # create a logical vector that indicates whether the current sequence name matches the species name.
  colname <- paste0(species) # create a column in the binary data frame for the current species name.
  DM_hm_nr_binary_df[colname] <- matches
}

# set the row names of the binary data frame to the sequence names.
rownames(DM_hm_nr_binary_df) <- DM_hm_nr_binary_df$seqname 

# remove the sequence name column from the binary data frame.
DM_hm_nr_binary_df <- DM_hm_nr_binary_df[,-c(1)]

# replace all NA values in the binary data frame with 0.
DM_hm_nr_binary_df[is.na(DM_hm_nr_binary_df)] <- 0

# optional step: transpose the binary data frame and convert it back to a data frame. this changes the place of rows and columns.
tDM_hm_nr_binary_df <- t(DM_hm_nr_binary_df) %>% as.data.frame()

```

##### Per Genome:
The code below also creates a binary data frame that indicates whether each sequence in the `DM_hm_nr_meta` and `FLS2_meta` data frames belongs to a particular organism (genome) in [`org_list`](trees/org_list.txt). The line-by-line annotation follows the previous script but based on organism names.

```R
# PERU
DM_hm_nr_org_binary_df <- data.frame(seqname = DM_hm_nr_meta$seqname)

for (i in 1:nrow(org_list)) {
  organism <- org_list$Organism[i]
  matches <- ifelse(DM_hm_nr_meta$Organism == organism, 1, 0)
  colname <- paste0(organism)
  DM_hm_nr_org_binary_df[colname] <- matches
}

rownames(DM_hm_nr_org_binary_df) <- DM_hm_nr_org_binary_df$seqname

DM_hm_nr_org_binary_df <- DM_hm_nr_org_binary_df[,-c(1)]

DM_hm_nr_org_binary_df[is.na(DM_hm_nr_org_binary_df)] <- 0


# FLS2
FLS2_org_binary_df <- data.frame(seqname = FLS2_meta$seqname)

for (i in 1:nrow(org_list)) {
  organism <- org_list$Organism[i]
  matches <- ifelse(FLS2_meta$Organism == organism, 1, 0)
  colname <- paste0(organism)
  FLS2_org_binary_df[colname] <- matches
}

rownames(FLS2_org_binary_df) <- FLS2_org_binary_df$seqname

FLS2_org_binary_df <- FLS2_org_binary_df[,-c(1)]

FLS2_org_binary_df[is.na(FLS2_org_binary_df)] <- 0

```

- Generating the plots
##### Figure 4A:
The `pheatmap` function is used to create heatmaps. The resulting heatmaps shows the distribution of the binary values in the `DM_hm_nr_binary_df` data frame. The rows of the heatmaps represent the species in the data frame and the columns represent the sequences. The colors in the heatmaps represent the presence (dark grey) or absence (white) in the data frame.

```R
library(pheatmap)

# custom color palette used for the plots
mypalette <- colorRampPalette(c("#F2F4F4","#424949"))


# PERU clade (fig 4A)
p_DM_hm_nr <- pheatmap(tDM_hm_nr_binary_df, # this transposed version of the binary dataframe allows to visualize the heatmap horizontally
                       cluster_rows = FALSE, cluster_cols = FALSE, 
                       show_colnames = TRUE, show_rownames = TRUE,
                       color = mypalette(2), 
                       border_color = "white", 
                       cellwidth = 10, cellheight = 10, # optional (as in Fig 4A)
                       gaps_row = (1:29)) # optional (as in Fig 4A)

ggsave(plot = p_DM_hm_nr, filename = "DM_hm_nr_pa.pdf",width = 35, height = 15, units = "in", dpi = "retina", device = "pdf", 
       path = "/path/to/your/directory/")
      
```

The following code snippets also create the heatmaps that appear in Figure S10, but instead of using species names, they use organism names.

##### Figure S10A:
```R
p_DM_hm_nr_org <- pheatmap(DM_hm_nr_org_binary_df, 
                        cluster_rows = FALSE, cluster_cols = FALSE, 
                        show_colnames = TRUE, show_rownames = TRUE,
                        color = mypalette(2), 
                        border_color = "white", 
                        cellwidth = 10, cellheight = 10) # optional (as in fig S10A)

# optional: you can export the plot
ggsave(plot = p_DM_hm_nr_org, filename = "DM_hm_nr_pa_org.pdf",width = 25, height = 35, units = "in", dpi = "retina", device = "pdf", path = "/path/to/your/directory/")

```

##### Figure S10B:

```R
p_FLS2_org <- pheatmap(FLS2_org_binary_df, 
                     cluster_rows = FALSE, cluster_cols = FALSE, 
                     show_colnames = TRUE, show_rownames = TRUE,
                     color = mypalette(2), 
                     border_color = "white", cellwidth = 10, cellheight = 10)

# optional: you can export the plot
ggsave(plot = p_FLS2_org, filename = "FLS2_pa_org.pdf",width = 25, height = 25, units = "in", dpi = "retina", device = "pdf", 
       path = "/path/to/your/directory/")

```
