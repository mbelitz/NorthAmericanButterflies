library(dplyr)
library(ggtree)
library(ggplot2)
# read in traits
st <- read.csv("MikeDir/data/TraitsButterflyCrossley.csv")
# read in phylogeny
# read in phylogeny
t <- read.tree(file = "MikeDir/data/iScience/SupDryad_treepl.tre")
length(t$tip.label)
head(t$tip.label)
# read in inat species
spp <- read.csv("MikeDir/data/SimpletonsSpeciesTraits.csv")
spp <- spp %>% 
  mutate(name = paste(Genus..iNat., species..iNat.))

spp <- spp %>% 
  mutate(name = if_else(condition =  name == " ", 
                        true = paste(NABA.Genus, NABA.species),
                        false = name))

# read in tiplabels that will match Spp of mdf
tip_labels_df <- read.csv("MikeDir/data/iScience/tipLabel_incidental.csv") # all tip labels
head(tip_labels_df$validName)

t$tip.label <- tip_labels_df$validName
head(tip_labels_df$validName)
head(t$tip.label)


tip_labels_df_lj <- left_join(tip_labels_df, spp, by = c("validName" = "name")) %>% 
  filter(!is.na(Simpleton.grouping.code))

# great, we have traits for all species so let's make into mdf
mdf <- tip_labels_df_lj %>% 
  filter(Simpleton.grouping.code != "OTH") %>% 
  filter(Simpleton.grouping.code != "Oth")

# drop stuff not in phylogy/analysis
mdf_phylo <- mdf %>% 
  rename(Species = validName) %>% 
  filter(Species %in% t$tip.label)

tree_sp <- t$tip.label
sppNotInAnalysis <- data.frame(Species = tree_sp) %>% 
  filter(!Species %in% mdf_phylo$Species)

t2 <- ape::drop.tip(t, tip = sppNotInAnalysis$Species)
t2 <- ape::drop.tip(t2, tip = "Lycaeides melissa")

p <- ggtree(t2, layout='circular')

phylo_df <- as.data.frame(p$data) 
phylo_df <- left_join(phylo_df, tip_labels_df, by = c("label" = "validName"))
phylo_df <- phylo_df%>% 
  filter(!is.na(tip_label))
library(stringr)
phylo_df <- phylo_df %>% 
  mutate(family = stringr::word(tip_label, 1,1,sep = fixed("_")))
phylo_df <- left_join(mdf_phylo, phylo_df, by = "tip_label")




p +
  geom_tile(data = phylo_df, aes(x = 135, y = y, fill = Simpleton.grouping.code),
            width = 20, color = "white",
            inherit.aes = F) + 
  geom_point(data = phylo_df, aes(x = x, y = y, color = family), size = 2) +
  scale_fill_viridis_d(option = "mako") +
  scale_color_viridis_d(option = "inferno") +
  labs(fill = "Overwintering stage", color = "Family")

ggsave(filename = "MikeDir/figures/phylogeny_incidentalDataSpecies.png", width = 14, height = 10)
