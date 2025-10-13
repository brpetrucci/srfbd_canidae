###########################################
##     Data cleanup for Canidae data     ##
## Chapter 2 - Bruno do Rosario Petrucci ##
###########################################
 
###
# load packages

# ape
library(ape)

# PBDB
library(paleobioDB)

# palaeoverse
library(palaeoverse)

# stringr
library(stringr)

###
# read data

# base directory
base_dir <- "/Users/petrucci/Documents/research/srfbd_canidae/empirical/"

# raw data directory
raw_dir <- paste0(base_dir, "raw_data/")

# fossil occurrences
occurrences <- read.csv(paste0(raw_dir, "occurrences_updated.csv"))[, -1]
occurrences <- occurrences[!is.na(occurrences$taxon), ]

# remove all occurrences with >15my resolution
occurrences <- occurrences[!(occurrences$early_age - occurrences$late_age > 15), ]

# molecular data
mol <- read.nexus.data(paste0(raw_dir, "mol.nex"))

# extant taxa
ext_taxa <- names(mol)[-which(names(mol) %in% c("Aenocyon_dirus", 
                                                "Dusicyon_australis"))]

# morphological data
morpho <- read.nexus.data(paste0(raw_dir, "morpho.nex"))

# function to check which characters have no variation
no_var_char <- function(morpho) {
  # make a results vector
  res <- c()
  
  for (i in 1:length(morpho[[1]])) {
    # get list of values
    vals <- unlist(lapply(morpho, function(x) x[i]))
    
    # get unique characters present
    chars <- unique(vals)
    
    # if there is only one non-? character, add to res
    if (sum(chars != "?") == 1) {
      res <- c(res, i)
    }
  }
  
  # return res
  return(res)
}


# function to switch 0/1 to {01}
rewrite_polymorphisms <- function(morpho) {
  # create new morpho list
  new_morpho <- vector("list", length(morpho))
  
  # iterate through morpho
  for (i in seq_along(morpho)) {
    # modify it
    new_morpho[[i]] <- ifelse(grepl("^(\\d+(?:/\\d+)+)$", 
                                    morpho[[i]], perl = TRUE),
                              paste0("{", gsub("/", "", morpho[[i]]), "}"),
                              morpho[[i]])
  }
  
  # name it
  names(new_morpho) <- names(morpho)
  
  # return new morpho
  return(new_morpho)
}

# get the characters with no variation
no_var <- no_var_char(morpho)

# remove the characters with no variation
morpho <- lapply(morpho, function(x) x[-no_var])

# rewrite polymorphisms
morpho <- rewrite_polymorphisms(morpho)

# get taxon names
taxa_names <- unique(c(names(morpho), names(mol)))

###
# extract range data

# function returning max and min for each range for a given taxon name
get_range <- function(occurrences, name, ext_as_zero = TRUE) {
  # get all occurrences with that name
  named_occs <- occurrences[occurrences$taxon == name, ]
  
  # number of occurrences
  n_occs <- nrow(named_occs)
  
  # get highest FA occurrences
  fa_high <- named_occs[which(named_occs$early_age ==
                                max(named_occs$early_age)), ]
  
  # max and min for fa
  fa_max <- max(fa_high$early_age)
  fa_min <- max(fa_high$late_age)
  
  # get lowest last age occurrence
  la_low <- named_occs[which(named_occs$early_age ==
                               min(named_occs$early_age)), ]
  
  # max and min for la
  la_min <- min(la_low$late_age)
  la_max <- max(la_low$early_age)
  
  return(c(name, fa_max, fa_min, la_max, la_min, n_occs))
}

# function to build ranges
build_ranges <- function(occurrences, taxa_names) {
  # apply to all species
  ranges <- lapply(unique(occurrences$taxon), 
                   function(x) get_range(occurrences, x))
  
  # make it a data frame
  ranges_df <- t(as.data.frame(ranges))
  colnames(ranges_df) <- c("taxon", "fa_max", "fa_min", "la_max", "la_min", "k")
  
  # add data without fossil occurrences to ranges
  nofossil_taxa <- taxa_names[!(taxa_names %in% ranges_df[, 1])]
  ranges_df <- rbind(ranges_df, data.frame(taxon = nofossil_taxa,
                                           fa_max = rep(0, length(nofossil_taxa)),
                                           fa_min = rep(0, length(nofossil_taxa)),
                                           la_max = rep(0, length(nofossil_taxa)),
                                           la_min = rep(0, length(nofossil_taxa)),
                                           k = 0))
  rownames(ranges_df) <- 1:nrow(ranges_df)
  
  # final ranges (change extants to 0 0)
  ranges_final <- ranges_df
  
  # iterate through extant taxa
  ranges_final$la_max[ranges_final$taxon %in% ext_taxa] <- 0
  ranges_final$la_min[ranges_final$taxon %in% ext_taxa] <- 0
  
  # reorder it alphabetically
  ranges_final <- ranges_final[order(ranges_final$taxon), ]
  rownames(ranges_final) <- 1:nrow(ranges_df)
  
  # return ranges_final
  return(ranges_final)
}

# get ranges for full dataset
ranges <- build_ranges(occurrences, taxa_names)

# write ranges to a file
write.table(ranges, paste0(base_dir, "data/ranges.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# add ?? for molecular data for species without molecular data
nomol_taxa <- taxa_names[!(taxa_names %in% names(mol))]
nomol_mol <- lapply(nomol_taxa, function(x) rep("?", length(mol[[1]])))
names(nomol_mol) <- nomol_taxa
mol_complete <- c(mol, nomol_mol)

# reorder mol alphabetically
mol_complete <- mol_complete[order(names(mol_complete))]

# write DNA to a file
write.nexus.data(mol_complete, paste0(base_dir, "data/mol.nex"))

# reading partitions and writing to DNA file
partitions <- readLines(paste0(raw_dir, "partitions_simplified.txt"))
write(partitions, paste0(base_dir, "data/mol.nex"), append = TRUE)

# add ?? for morphological data for species without morphological data
nomorpho_taxa <- taxa_names[!(taxa_names %in% names(morpho))]
nomorpho_morpho <- lapply(nomorpho_taxa, function(x) rep("?", length(morpho[[1]])))
names(nomorpho_morpho) <- nomorpho_taxa
morpho_complete <- c(morpho, nomorpho_morpho)

# reorder morpho alphahetically
morpho_complete <- morpho_complete[order(names(morpho_complete))]

# write morpho to a file
write.nexus.data(morpho_complete, paste0(base_dir, "data/morpho.nex"),
                 format = "standard")

###
# function to save subset of data
subset_data <- function(mol, morpho, taxa) {
  # get subset mol and morpho data
  mol_subset <- mol[taxa]
  morpho_subset <- morpho[taxa]
  
  # find which morphological characters are now constant
  no_var_subset <- no_var_char(morpho_subset)
  
  # select only the characters that vary
  if (length(no_var_subset) > 0)
    morpho_subset <- lapply(morpho_subset, function(x) x[-no_var_subset])
  
  # return mol and morpho for the subset
  return(list(MOL = mol_subset, MORPHO = morpho_subset))
}

###
# creating a data set without occurrences with >5my uncertainty

# filter occurrences
occurrences_filter <- (occurrences$early_age - occurrences$late_age) > 5
ftrd_occs <- occurrences[!occurrences_filter, ]

# taxa names for ftrd occurrences
removed_taxa <- unique(occurrences$taxon[!(occurrences$taxon %in% ftrd_occs$taxon)])
ftrd_taxa_names <- unique(names(mol_complete)[!(names(mol_complete) %in% removed_taxa)])

# get data for ftrd taxa
ftrd_data <- subset_data(mol_complete, morpho_complete, ftrd_taxa_names)
ftrd_morpho <- ftrd_data$MORPHO
ftrd_mol <- ftrd_data$MOL

# filter ranges
ftrd_ranges <- build_ranges(ftrd_occs, names(ftrd_mol))

# write ranges and morpho data
write.table(ftrd_ranges, paste0(base_dir, "data/ftrd_ranges.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)
write.nexus.data(ftrd_morpho, paste0(base_dir, "data/ftrd_morpho.nex"),
                 format = "standard")
write.nexus.data(ftrd_mol, paste0(base_dir, "data/ftrd_mol.nex"))
write(partitions, paste0(base_dir, "data/ftrd_mol.nex"), append = TRUE)

###
# saving extant and recent data

# get data for extant taxa
extant_data <- subset_data(mol_complete, morpho_complete, ext_taxa)
extant_morpho <- extant_data$MORPHO
extant_mol <- extant_data$MOL

# and for recent taxa
recent_data <- subset_data(mol_complete, morpho_complete, 
                           c(ext_taxa, "Aenocyon_dirus", "Dusicyon_australis"))
recent_morpho <- recent_data$MORPHO
recent_mol <- recent_data$MOL

# save extant data
write.nexus.data(extant_morpho, paste0(base_dir, "data/ext_morpho.nex"),
                 format = "standard")
write.nexus.data(extant_mol, paste0(base_dir, "data/ext_mol.nex"))
write(partitions, paste0(base_dir, "data/ext_mol.nex"), append = TRUE)

# and recent
write.nexus.data(recent_morpho, paste0(base_dir, "data/recent_morpho.nex"),
                 format = "standard")
write.nexus.data(recent_mol, paste0(base_dir, "data/recent_mol.nex"))
write(partitions, paste0(base_dir, "data/recent_mol.nex"), append = TRUE)

###
# plot ranges

# organize it
ranges_plot_df <- data.frame(x1 = as.numeric(ranges$fa_max), 
                             x2 = as.numeric(ranges$la_max),
                             y = 1:nrow(ranges))

# plot
plot(1, type = "n", xlab = "", axes = FALSE,
     ylab = "", xlim = c(40, 0),  
     ylim = c(0, 158)) 
segments(x0 = ranges_plot_df$x1, y0 = ranges_plot_df$y, 
         x1 = ranges_plot_df$x2, y1 = ranges_plot_df$y,
         lwd = 2)
axis_geo(side = 1, intervals = "epochs")
