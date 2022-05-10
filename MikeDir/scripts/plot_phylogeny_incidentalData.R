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
  mutate(name = paste(Genus..iNat., species..iNat.)) %>% 
  filter(name != " ") %>% 
  select(name, Simpleton.grouping.code)

# read in tiplabels that will match Spp of mdf
tip_labels_df <- read.csv("MikeDir/data/iScience/tipLabel.csv") # all tip labels
head(tip_labels_df$validName)

t$tip.label <- tip_labels_df$validName
head(tip_labels_df$validName)
head(t$tip.label)


tip_labels_df_lj <- left_join(tip_labels_df, spp, by = c("validName" = "name")) %>% 
  filter(!is.na(Simpleton.grouping.code))

# great, we have traits for all species so let's make into mdf
mdf <- tip_labels_df_lj %>% 
  filter(Simpleton.grouping.code != "OTH")

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

phylo_df <- data.frame(species = t2$tip.label, node = t2$Nnode)
t3 <- read.csv("MikeDir/data/t3.csv")

p +
  geom_tile(data = t3, aes(x = 135, y = y, fill = Simpleton.grouping.code),
            width = 20, color = "white",
            inherit.aes = F) + 
  geom_point(data = t3, aes(x = x, y = y, color = family), size = 2) +
  scale_fill_viridis_d(option = "mako") +
  scale_color_viridis_d(option = "inferno") +
  labs(fill = "Overwintering stage", color = "Family")

ggsave(filename = "MikeDir/figures/phylogeny_incidentalDataSpecies.png", width = 14, height = 10)
