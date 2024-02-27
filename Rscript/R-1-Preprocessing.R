################################################
### Preprocessing follows
### https://github.com/chiras/metabarcoding_pipeline
### https://doi.org/10.1098/rstb.2021.0171

print("Data loading and preparation")

################################################
# Loading in data
data.tax <- tax_table(as.matrix(read.table("Data/total.taxonomy.sed.vsearch", header=T,row.names=1,fill=T,sep=",")))

## Community table
data.otu <- otu_table(read.table("Data/total.asv_table.bytable.sed.txt"), taxa_are_rows=T)

## Sample metadata and adapt names
data.map <- sample_data(read.table("Data/Map4.txt", header=T, row.names=1,  sep="\t", fill=T))
sample_names(data.map) <- gsub("-",".",sample_names(data.map))

## check metadata vs. samples in sequencing data consistency
sample_names(data.map )[!(sample_names(data.map ) %in% sample_names(data.otu))]
sample_names(data.otu )[!(sample_names(data.otu ) %in% sample_names(data.map))]

## merge the three tables to a single phylseq object
(data.comp <- merge_phyloseq(data.otu,data.tax,data.map ))

## information about crop status
crops = read.table("Data/crops.list")
crops$V1 = gsub("_"," ",crops$V1)

################################################
# Preprocessing
## given hierarchical classification options at the end, we have to propagate the taxonomy over taxonomic levels to not throw out sequences only classified to higher tax levels
data.comp <- propagate_incomplete_taxonomy(data.comp)

## filtering irrelevant taxa, zb. unresolved, algae, fungi etc
data.comp.filter <- remove_unresolved_taxa(data.comp)

## Make taxa labels nice for plots
data.comp.filter <- replace_tax_prefixes(data.comp.filter)

## Multiple ASVs might represent the same species, here they are collated
(data.species <- tax_glom(data.comp.filter,taxrank="species"))
taxa_names(data.species) <- tax_table(data.species)[,"species"]

## Transform to relative data
data.species.rel = transform_sample_counts(data.species, function(x) x/sum(x))

## low abundance filtering
otu_table(data.species.rel)[otu_table(data.species.rel)<0.01 ]<-0
otu_table(data.species)[otu_table(data.species.rel)<0.01 ]<-0
data.species.filter		= prune_taxa(taxa_sums(data.species)>0, data.species)
(data.species.rel.filter = prune_taxa(rowSums(otu_table(data.species.rel))>0, data.species.rel))

## define parameters for plot definitions
ntaxa <- length(taxa_names(data.species.rel.filter))
