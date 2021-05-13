library(hash)
library(dplyr)

# Function which reads in the Blosum txt file to creating a scoring dictionary
processFile = function(blos_path) {
  con = file(blos_path, "r")
  
  line = readLines(con, n = 1)
  letters = strsplit(line, "\\s+")[[1]]
  
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    
    line = strsplit(line, "\\s+")[[1]]
    key = line[[1]]
    temp_dict = hash() 
    for (idx in 2:24) {
      temp_dict[[letters[idx]]] <- strtoi(line[[idx]]) - 11
    }
    scoreDict[[key]] = temp_dict
  }
  close(con)
}

# Function which uses the scoring dictionary to provide the A/A changes score
getScore <- function(x) {
  a = substring(x, 1, 1) 
  b = substring(x, 3, 3)
  return(scoreDict[[a]][[b]])
}

# Function which converts the verbal Impact scores into numbers
getImpactScore <- function(x) {
  pos = c("MODIFIER", "LOW", "MODERATE", "HIGH")
  idx = which(x == pos)[[1]]
  return(idx - 1)
}


#CODE TO READ IN TXT OUTPUT FROM VEP which gives us an annotated version of the uploaded file with all gene names
df = read.delim("./data/1AkM6SVSzu8HTWfp.txt", header = TRUE, stringsAsFactors = FALSE, quote = "")

# Taking a peek at the data table
head(df)
# What are the possible values for consequence? 
unique(df$Consequence)

# Filter out undesired mutations to get top genes, but we will want to reference original 
# df later to get all mutations
desired = list("missense_variant", "missense_variant,NMD_transcript_variant",
             "missense_variant,splice_region_variant,NMD_transcript_variant",
             "start_lost", "start_lost,NMD_transcript_variant",
             "missense_variant,splice_region_variant", "stop_retained_variant",
             "stop_lost", "stop_gained")

df_filtered = subset(df, Consequence %in% desired)
head(df_filtered)

# Find 10 highest appearing gene names in filtered table
sort(table(df_filtered$Gene), decreasing = TRUE)[1:10]

# Add a score based on the AA changes
blos_path = "./data/blosum62.txt"
scoreDict = hash() 

processFile(blos_path)

# Add column for scores
df_filtered$Score = unlist(lapply(df_filtered$Amino_acids, function(x){ifelse(is.na(x),NA, getScore(x))}))
head(df_filtered)

# Wrangle sorting so that we get the lowest scoring gene
df_sum <-  
  group_by(df_filtered, Gene) %>% 
  summarise(Score = sum(Score))

# Going to grab the mutations for gene ENSG00000051180
filter(df, Gene=='ENSG00000051180')

# Grab the original VEP run data
df2 = read.delim("./data/RPwvCHag6FBO2JRv.txt", header = TRUE, stringsAsFactors = FALSE, quote = "")
head(df2)
filter(df2, Gene=='ENSG00000051180' & Codons != '-')

df2_group = group_by(df2, Gene)
head(df2_group)
glimpse(df2_group)

# Try to get unique protein positions per gene that are mutated
df2_groups =
  group_by(df2, Gene) %>%
  summarise(distinct_visit_ids = n_distinct(Protein_position))

arrange(df2_groups, desc(distinct_visit_ids))

# Find 10 highest appearing gene names in filtered table
sort(table(df2$Gene), decreasing = TRUE)[1:10]
sort(table(filter(df2, IMPACT == 'HIGH')$Gene), decreasing = TRUE)

df2$ImpactScore = unlist(lapply(df2$IMPACT, function(x){ifelse(is.na(x), NA, getImpactScore(x))}))

df2_groups =
  group_by(df2, Gene) %>%
  summarise(TotalImpactScore = sum(ImpactScore))
arrange(df2_groups, desc(TotalImpactScore))

df_sum =
  group_by(df_filtered, Gene) %>%
  summarise(TotalScore = sum(unlist(Score)))
arrange(df_sum, TotalScore)

# Take a look at AKAP13 to get all protein mutations
filter(df2, Gene=='ENSG00000170776' & Codons != '-' & Consequence != 'synonymous_variant')
