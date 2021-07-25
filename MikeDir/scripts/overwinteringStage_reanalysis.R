#load libraries
library(lme4)
library(lmerTest)
library(car)

bfly_trends = read.table("data/butterfly_traits_envars_trends_50km_NoMergeJul_m5_trim_6traits.txt",sep='\t',as.is=T,check.names=F,header=T)

# build first model

m <- lmer(Abundance.trend ~ Cropland.trend + Built.2005_2015 + Precip.1993_2018 +
          Temp.1993_2018 + LarvalColor + LarvalHair + AdultSize + AdultColor +  
            Diet.breadth.families + (1 | grid_id) + (1 | Species),
          data = bfly_trends, REML = F, lmerControl(optimizer = "bobyqa"))

ms <- step(m)

## so our results are similar to crossley et all, which is cool 
m_top <- lmer(Abundance.trend ~ Precip.1993_2018 + Temp.1993_2018 + 
                (1 | grid_id) + (1 | Species),
              data = bfly_trends, REML = F, lmerControl(optimizer = "bobyqa"))

summary(m_top)

## let's get a list of species
spp_list <- bfly_trends$Species %>% unique()

write.csv(spp_list, "MikeDir/data/spp_list.csv", row.names = F)

#simpleton traits
st <- read.csv("MikeDir/data/SimpletonsSpeciesTraits.csv")
st <- st %>% 
  mutate(Species = paste(NABA.Genus,NABA.species))


bfly_trends_st <- left_join(bfly_trends, st)

mdf <- bfly_trends_st %>% 
  filter(!is.na(Simpleton.grouping.code)) %>% 
  mutate(Simpleton.grouping.code = 
           if_else(condition = Simpleton.grouping.code == "Oth", 
                 true = "OTH", false = Simpleton.grouping.code))

m <- lmer(Abundance.trend ~ Cropland.trend + Built.2005_2015 + Precip.1993_2018 +
            Temp.1993_2018 + LarvalColor + LarvalHair + AdultSize + AdultColor +  
            Diet.breadth.families + Simpleton.grouping.code +
            Temp.1993_2018:Simpleton.grouping.code +
            (1 | grid_id) + (1 | Species),
          data = mdf, REML = F, lmerControl(optimizer = "bobyqa"))

ms <- step(m)

m_top <- lmer(Abundance.trend ~ Precip.1993_2018 + Temp.1993_2018 + Simpleton.grouping.code + 
                Temp.1993_2018:Simpleton.grouping.code +
                (1 | grid_id) + (1 | Species),
              data = mdf, REML = F, lmerControl(optimizer = "bobyqa"))

summary(m_top)                
vif(m_top)                

library(sjPlot)
plot_model(m_top, type = "pred", terms = "Simpleton.grouping.code")
plot_model(m_top, type = "pred", terms = c("Temp.1993_2018", "Simpleton.grouping.code"))
