# larsen et al. extent

#load libraries
library(lme4)
library(lmerTest)
library(car)
library(tidyverse)
library(MuMIn)
library(ape)
library(stringr)
library(sjPlot)

# read in crossley abundance trend data
bfly_trends = read.table("data/butterfly_traits_envars_trends_50km_NoMergeJul_m5_trim_6traits.txt",
                         sep='\t',as.is=T,check.names=F,header=T)

# first we shall look at butterflies across north america results with traits we think are important
#read in traits
st <- read.csv("MikeDir/data/TraitsButterflyCrossley.csv")
# join to butterfly trend data
bfly_trends_st <- left_join(bfly_trends, st)

# do any species not have OWS?
dplyr::filter(bfly_trends_st, is.na(overwintering.stage))
# voltinism?
dplyr::filter(bfly_trends_st, is.na(voltinism))
# we have traits for every species
# what sort of unique traits do we have?
unique(bfly_trends_st$voltinism) #F for facultative, "M" for multivoltine, and "U" for univoltine "" needs to be removed
unique(bfly_trends_st$overwintering.stage) # so lots of stuff, going to make "A (migrate)" into "M" and then 

# great, we have traits for all species so let's make into mdf
mdf <- bfly_trends_st %>% 
  mutate(overwintering.stage = if_else(condition = overwintering.stage == "A (migrate)",
                                       true = "M",
                                       false = overwintering.stage)) %>% 
  filter(overwintering.stage %in% c("L", "A", "P", "M", "E")) %>% 
  filter(voltinism != "")

# let's filter to our study extent
# (32° N, 50° N) and (95° W, 60° W), 
mdf <- mdf %>% 
  filter(Lon > -95 & Lon < -60 ) %>% 
  filter(Lat > 32 & Lat < 50)
# filter to only species that overwinter as eggs, larvae and pupae
mdf <- mdf %>% 
  filter(overwintering.stage %in% c("E", "L", "P"))

# and now filter only to species in our modeling framework
naba_larsen <- read.csv("MikeDir/data/NABA_abund_taxa.csv")
#join larsen dataset to mdf
crossley_spp <- data.frame(Species = unique(mdf$Species))
lj <- left_join(crossley_spp, naba_larsen, by = c("Species" = "ScientificName"))
# so i think everything missing is not in larsen dataset not a taxonomy problem
mdf2 <- left_join(mdf, naba_larsen, by = c("Species" = "ScientificName")) %>% 
  filter(!is.na(mean.abund))

# down to 4599 rows
# modeling time
m <- lmer(Abundance.trend ~ Cropland.trend + Built.2005_2015 + Precip.1993_2018 +
            Temp.1993_2018 + LarvalColor + LarvalHair + AdultSize + AdultColor +  
            Diet.breadth.families + overwintering.stage + voltinism + 
            Temp.1993_2018:overwintering.stage + 
            Precip.1993_2018:overwintering.stage +
            (1 | grid_id) + (1 | Species),
          data = mdf, REML = F, lmerControl(optimizer = "bobyqa"))

ms <- step(m)
ms

m_top <- lmer(Abundance.trend ~ Precip.1993_2018 + Temp.1993_2018 + 
                overwintering.stage + 
                Precip.1993_2018:overwintering.stage +
                (1 | grid_id) + (1 | Species),
              data = mdf, REML = F, lmerControl(optimizer = "bobyqa"))


summary(m_top)                
vif(m_top)                
r.squaredGLMM(m_top) ## 0.0368 # 0.6301

# plot results
plot_model(m_top, type = "pred", terms = "overwintering.stage")
plot_model(m_top, type = "pred", terms = c("Precip.1993_2018", "overwintering.stage"))

#
## now take that top model and make it a phylogenetic linear mixed model
# read in phylogeny
t <- read.tree(file = "MikeDir/data/iScience/SupDryad_treepl.tre")
length(t$tip.label)
head(t$tip.label)
# read in tiplabels that will match Spp of mdf
tip_labels_df <- read.csv("MikeDir/data/iScience/tipLabel.csv") # all tip labels
head(tip_labels_df$validName)

t$tip.label <- tip_labels_df$validName
head(tip_labels_df$validName)
head(t$tip.label)

# drop species not in phylogeny
mdf_phylo <- mdf2 %>% 
  filter(Species %in% t$tip.label)

# now we drop our data from 9535 rows of Crossley to 4471 
# total of 52.6% loss

## now add in a phylogenetic effect to the modeling framework
library(phyr)

mdf_phylo <- mdf_phylo %>% 
  mutate(Precip.1993_2018 = scale(Precip.1993_2018),
         Temp.1993_2018 = scale(Temp.1993_2018))

tree_sp <- t$tip.label
sppNotInAnalysis <- data.frame(Species = tree_sp) %>% 
  filter(!Species %in% mdf_phylo$Species)

t2 <- ape::drop.tip(t, tip = sppNotInAnalysis$Species)

pglmm <- pglmm(Abundance.trend ~ Precip.1993_2018 + Temp.1993_2018 + 
                 overwintering.stage + 
                 Precip.1993_2018:overwintering.stage +
                 (1 | grid_id) + (1 | Species__) ,
               data = mdf_phylo,
               cov_ranef = list(Species = t2), 
               bayes = T, verbose = T)

# look into model summary
summary(pglmm)
# look at posterior distributions
p <- phyr::plot_bayes(x = pglmm)

ggsave(plot = p, filename = "MikeDir/figures/pglmm_effects_LarsenExtent.png",
       dpi = 400, width = 12, height = 8)

x = pglmm
n_samp = 1000
sort = TRUE
re.names <- names(x$random.effects)
if (x$family == "gaussian") re.names <- c("residual", re.names)
random_samps <- lapply(x$inla.model$marginals.hyperpar, 
                       function(x) INLA::inla.rmarginal(n_samp, INLA::inla.tmarginal(function(x) sqrt(1 / x), x))) %>%
  setNames(re.names) %>%
  dplyr::as_tibble() %>%
  tidyr::pivot_longer(cols = dplyr::everything(),
                      names_to = "var",
                      values_to = "val") %>%
  dplyr::mutate(effect_type = "Random Effects")

fixed_samps <- lapply(x$inla.model$marginals.fixed, function(x) INLA::inla.rmarginal(n_samp, x)) %>%
  dplyr::as_tibble() %>%
  tidyr::pivot_longer(cols = dplyr::everything(),
                      names_to = "var",
                      values_to = "val") %>%
  dplyr::mutate(effect_type = "Fixed Effects")

samps <- dplyr::bind_rows(random_samps, fixed_samps) %>%
  dplyr::mutate(effect_type = factor(effect_type, 
                                     levels = c("Random Effects", "Fixed Effects")))

ci <- samps %>%
  dplyr::group_by(var, effect_type) %>%
  dplyr::summarise(lower = quantile(val, 0.025),
                   upper = quantile(val, 0.975),
                   mean = mean(val),
                   .groups = "drop_last")

if(sort){
  ci <- dplyr::arrange(ci, mean) %>% dplyr::ungroup() %>% 
    dplyr::mutate(var = factor(as.character(var), levels = as.character(var)))
}

sig_vars <- ci %>%
  dplyr::mutate(sig = ifelse(effect_type == "Random Effects",
                             "CI no overlap with zero",
                             ifelse(sign(lower) == sign(upper),
                                    "CI no overlap with zero",
                                    "CI overlaps zero"))) %>%
  dplyr::select(var, sig)

if(sort){
  samps <- dplyr::mutate(samps, var = factor(var, levels = levels(sig_vars$var)))
}

samps <- samps %>%
  dplyr::left_join(sig_vars, by = "var") %>%
  dplyr::group_by(var) %>%
  dplyr::filter(abs(val - mean(val)) < (10 * sd(val))) %>% 
  dplyr::ungroup()

pal <- c("#fc8d62", "#8da0cb")
p <- ggplot2::ggplot(filter(samps, effect_type == "Random Effects"),
                     ggplot2::aes(val, var, height = ..density..)) +
  ggridges::geom_density_ridges(ggplot2::aes(alpha = sig, fill = sig), 
                                stat = "density", adjust = 2, color = "gray70") +
  ggplot2::geom_point(ggplot2::aes(x = mean, y = var),
                      data = filter(ci, effect_type == "Random Effects"), inherit.aes = FALSE) +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower, xmax = upper, y = var),
                          data = filter(ci, effect_type == "Random Effects"),
                          inherit.aes = FALSE, height = 0.1) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  ggplot2::scale_alpha_manual(values = c(0.8, 0.2)) +
  ggplot2::scale_fill_manual(values = rev(pal)) +
  ggplot2::ylab("") +
  ggplot2::xlab("Estimate") +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "none",
                 axis.text = ggplot2::element_text(size = 14),
                 strip.text = ggplot2::element_text(size = 16))

p


ggsave(plot = p, filename = "MikeDir/figures/pglmm_RandomEffects_LarsenExtent.png",
       dpi = 400, width = 12, height = 8)

# now plot interactions
pred_vals <- -3:3

inla_mdf <- mdf_phylo %>% 
  dplyr::select(Abundance.trend, Temp.1993_2018, 
                Precip.1993_2018, overwintering.stage,
                Species, grid_id)

term1 <- data.frame(
  Abundance.trend = rep(NA, length(pred_vals)),
  Precip.1993_2018 = pred_vals,
  Temp.1993_2018 = rep(NA, length(pred_vals)),
  overwintering.stage = rep("E", length(pred_vals)),
  grid_id = rep(NA, length(pred_vals)),
  Species = rep(NA, length(pred_vals)))

term2 <- data.frame(
  Abundance.trend = rep(NA, length(pred_vals)),
  Precip.1993_2018 = pred_vals,
  Temp.1993_2018 = rep(NA, length(pred_vals)),
  overwintering.stage = rep("L", length(pred_vals)),
  grid_id = rep(NA, length(pred_vals)),
  Species = rep(NA, length(pred_vals)))

term3 <- data.frame(
  Abundance.trend = rep(NA, length(pred_vals)),
  Precip.1993_2018 = pred_vals,
  Temp.1993_2018 = rep(NA, length(pred_vals)),
  overwintering.stage = rep("P", length(pred_vals)),
  grid_id = rep(NA, length(pred_vals)),
  Species = rep(NA, length(pred_vals)))


pred.df.1 <- rbind(inla_mdf, term1)
pred.df.2 <- rbind(inla_mdf, term2)
pred.df.3 <- rbind(inla_mdf, term3)

pglmm1 <- pglmm(Abundance.trend ~ Precip.1993_2018 + Temp.1993_2018 + 
                  overwintering.stage + 
                  Temp.1993_2018:overwintering.stage +
                  Precip.1993_2018:overwintering.stage +
                  (1 | grid_id) + (1 | Species__) ,
                data = pred.df.1,
                cov_ranef = list(Species = t2), 
                bayes = T, verbose = T)

pglmm2 <- pglmm(Abundance.trend ~ Precip.1993_2018 + Temp.1993_2018 + 
                  overwintering.stage + 
                  Temp.1993_2018:overwintering.stage +
                  Precip.1993_2018:overwintering.stage +
                  (1 | grid_id) + (1 | Species__) ,
                data = pred.df.2,
                cov_ranef = list(Species = t2), 
                bayes = T, verbose = T)

pglmm3 <- pglmm(Abundance.trend ~ Precip.1993_2018 + Temp.1993_2018 + 
                  overwintering.stage + 
                  Temp.1993_2018:overwintering.stage +
                  Precip.1993_2018:overwintering.stage +
                  (1 | grid_id) + (1 | Species__) ,
                data = pred.df.3,
                cov_ranef = list(Species = t2), 
                bayes = T, verbose = T)

## make df of results
rdf1 <- pglmm1$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Egg")

rdf2 <- pglmm2$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Larvae")

rdf3 <- pglmm3$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Pupae")

rdf_total <- rbind(rdf1, rdf2, rdf3)

prec_ows <- ggplot(rdf_total, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                            fill = V2), 
              alpha = 0.15) +
  geom_line(mapping = aes(color = V2), size = 1.05) +
  labs(x = "Precipitation", y = "Abundance trend", 
       fill = "Overwintering stage", color = "Overwintering stage") + 
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  theme_classic() +
  theme(legend.position = "bottom")

prec_ows

ggsave(plot = prec_ows, filename = "MikeDir/figures/prec_OWS_interaction_LarsenExtent.png",
       dpi = 450, width = 6, height = 3.5)

cp <- cowplot::plot_grid(prec_ows, p, labels = c("A", "B"))


ggsave(plot = cp, filename = "MikeDir/figures/LarsenExtent_Interactions&RandomEffect.png",
       dpi = 450, width = 9, height = 4)
