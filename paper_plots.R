library(ggplot2)
library(tidyr)
library(dplyr)
library(ggsurveillance)
library(scales)
library(egg)
library(ggh4x)
library(cowplot)

hex <- hue_pal()(3)

##### No cost Drug effects on death (not growth) rate, additive (not multiplicative) on MIC #####

# Function to calculate growth rate
growth_rate_genotype <- function(
    phi_max, phi_min, concentrations, MIC, 
    k, cost_vector, benefit_vector, genotype) {
  
  MIC_gen <- MIC + sum(benefit_vector*genotype)
  phi_max_gen <- phi_max - sum(cost_vector*genotype)
  
  phi <- phi_max_gen - (phi_max_gen - phi_min) * (concentrations / MIC_gen)^k / (
    (concentrations / MIC_gen)^k - phi_min / phi_max_gen
  )
  
  return(phi)
}


# Base parameters
phi_max <- 1
phi_min <- -1
MIC <- 20
k <- 3
concentrations <- 10^seq(-1, 3, length.out = 100)

low_ben = 80
high_ben = 200

#### Explain dose response curve ######

k=3
benefits_vector = c(low_ben, high_ben)
cost_vector <- c(0.2, 0.2)

all_genotypes = expand.grid(rep(list(c(0, 1)), 2))

explainDR_curve <- data.frame(
  rate = growth_rate_genotype(
    phi_max, phi_min, concentrations, MIC, 
    k, cost_vector, c(0,0), c(0,0)),
  genotype = '0 0', num_muts = 0) 

for (i in 2:nrow(all_genotypes)){
  geno = all_genotypes[i,]
  explainDR_curve = rbind(
    explainDR_curve, 
    data.frame(
      rate = growth_rate_genotype(
        phi_max, phi_min, concentrations, MIC, 
        k, cost_vector, benefits_vector, geno),
      genotype = paste(geno, collapse = ' '), num_muts = sum(geno!=0)
    )
  )
}

explainDR_curve$concentration = rep(
  concentrations, nrow(explainDR_curve)/length(concentrations))

explainDR_curve = explainDR_curve %>%
  dplyr::mutate(
    xmin_single = 10^(1.3),
    xmax_MIC_single = 10^(2.1),
    vlines = 10
    )

hex <- palette.colors(palette = "Okabe-Ito") #hue_pal()(7)

fix_spacing = 'genotype  '
# Create plot
ggplot(explainDR_curve, 
       aes(x=concentration, y=rate, color = genotype, 
           linetype = genotype, size=genotype)) +
  geom_rect(
    data = unique(explainDR_curve[c("xmin_single","xmax_MIC_single")]),
    aes(xmin =xmin_single, xmax = xmax_MIC_single,
        ymin=-Inf, ymax=Inf), 
    inherit.aes = FALSE, alpha=0.4, fill='darkgrey') +
  geom_rect(
    data = unique(explainDR_curve[c("xmin_single","vlines")]),
    aes(xmin =vlines, xmax = xmin_single,
        ymin=-Inf, ymax=Inf), 
    inherit.aes = FALSE, alpha=0.4, fill='orange') +
  # geom_vline(aes(xintercept = vlines), color = 'darkgrey', 
  #            linetype = 'dashed', linewidth = 0.8)+
  geom_line()+
  theme_bw(base_size = 23) +
  labs(color = 'genotype')+
  scale_color_manual(values = c('black',hex[2:5]))+
  scale_size_manual(values = c(0.7, 1.5, 1.5, 1.5))+
  scale_linetype_manual(values=c("solid", "solid", "solid", "dashed")) +
  ylab('Growth rate (arbitrary units)') +
  xlab('Drug concentration (arbitrary units)')  +
  scale_x_log10(label=label_power10(magnitude_only = TRUE), limits = c(1, 1e3))+
  scale_y_continuous(breaks = c(-1,0,1), minor_breaks = c(-0.75, -0.5, -0.25, 0.25, 0.5, 0.75)) +
  guides(linetype = guide_legend(override.aes = list(size = 4), keywidth = 4, keyheight = 1),
         colour=guide_legend(keywidth = 4, keyheight = 1))+
  labs(color = fix_spacing, linetype = fix_spacing, size = fix_spacing)

ggsave('/Users/louagema/Documents/Figures/CollaterEffectsConcept/explain_DR_curve.png',
       dpi = 300, height = 5, width = 7.5)


### Explain DR curve with CS ####

benefits_vector = c(low_ben, high_ben, -MIC/2)
cost_vector <- c(0.2, 0.2, 0.2)

all_genotypes = expand.grid(rep(list(c(0, 1)), 3))

explainDR_curve <- data.frame(
  rate = growth_rate_genotype(
    phi_max, phi_min, concentrations, MIC, 
    k, cost_vector, c(0,0,0), c(0,0,0)),
  genotype = '0 0 0 ', num_muts = 0) 

for (i in c(2,3,4,5)){
  print(i)
  geno = all_genotypes[i,]
  print(geno)
  explainDR_curve = rbind(
    explainDR_curve, 
    data.frame(
      rate = growth_rate_genotype(
        phi_max, phi_min, concentrations, MIC, 
        k, cost_vector, benefits_vector, geno),
      genotype = paste(geno, collapse = ' '), num_muts = sum(geno!=0)
    )
  )
}

explainDR_curve$concentration = rep(
  concentrations, nrow(explainDR_curve)/length(concentrations))

explainDR_curve = explainDR_curve %>%
  dplyr::mutate(
    xmin_single = 10^(1.3),
    xmax_MIC_single = 10^(2.1),
    vlines = 10
  )

hex <- palette.colors(palette = "Okabe-Ito") #hue_pal()(7)

fix_spacing = 'genotype  '
# Create plot
ggplot(explainDR_curve, 
       aes(x=concentration, y=rate, color = genotype, 
           linetype = genotype, size=genotype)) +
  geom_rect(
    data = unique(explainDR_curve[c("xmin_single","xmax_MIC_single")]),
    aes(xmin =xmin_single, xmax = xmax_MIC_single,
        ymin=-Inf, ymax=Inf), 
    inherit.aes = FALSE, alpha=0.4, fill='darkgrey') +
  geom_rect(
    data = unique(explainDR_curve[c("xmin_single","vlines")]),
    aes(xmin =vlines, xmax = xmin_single,
        ymin=-Inf, ymax=Inf), 
    inherit.aes = FALSE, alpha=0.4, fill='orange') +
  # geom_vline(aes(xintercept = vlines), color = 'darkgrey', 
  #            linetype = 'dashed', linewidth = 0.8)+
  geom_line()+
  theme_bw(base_size = 23) +
  labs(color = 'genotype')+
  scale_color_manual(values = c('black',hex[c(2,3,4,6)]))+
  scale_size_manual(values = c(0.7, 1.5, 1.5, 1.5, 1.5))+
  scale_linetype_manual(values=c("solid", "solid", "solid", "solid", "dashed")) +
  ylab('Growth rate (arbitrary units)') +
  xlab('Drug concentration (arbitrary units)')  +
  scale_x_log10(label=label_power10(magnitude_only = TRUE), limits = c(1, 1e3))+
  scale_y_continuous(breaks = c(-1,0,1), minor_breaks = c(-0.75, -0.5, -0.25, 0.25, 0.5, 0.75)) +
  guides(linetype = guide_legend(override.aes = list(size = 4), keywidth = 4, keyheight = 1),
         colour=guide_legend(keywidth = 4, keyheight = 1))+
  labs(color = fix_spacing, linetype = fix_spacing, size = fix_spacing)

ggsave('/Users/louagema/Documents/Figures/CollaterEffectsConcept/explain_DR_curveSandR.png',
       dpi = 300, height = 5, width = 7.5)

#### Overview figure #####

# Three loci with cost
cost_vector <- c(0, 0, 0)

# Same benefits, two loci
double_same <- c(high_ben, high_ben, 0)
single_same <- c(high_ben, 0, 0)

# One locus R
single_only <- c(high_ben, 0, 0)

# Low, high benefit, two loci
double_low_high <- c(low_ben, high_ben, 0)
single_low <- c(low_ben, 0, 0)
single_high <- c(0, high_ben, 0)

wildtype <- growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, c(0,0,0), c(0,0,0))

# Calculate growth rates for different scenarios
same_benefits_two <- growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, double_same, c(1,1,0))
same_benefits_one <- growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, single_same, c(1,0,0))

only_one <- growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, single_only, c(1,0,0))

low_high_two <- growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, double_low_high, c(1,1,0))
low_one <- growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, single_low, c(1,0,0))
high_one <- growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, single_high, c(0,1,0))

# Create dataframe (Note: you may need to adjust the variable names)
# The Python code references variables not defined in the notebook

names = c('same/similar benefits \n for beneficial mutations', 
          'only one beneficial \n mutation available',
          'low or high benefit \n for beneficial mutations')

nocost_death_additive <- data.frame(
  rate = rep(wildtype, 3),
  genotype = '[0, 0]', num_muts = '0',
  scenario = rep(names, each = length(wildtype))
) %>%
  rbind(data.frame(
    rate = same_benefits_one,
    genotype = '[1, 0]', num_muts = '1',
    scenario = names[1]
  )) %>%
  rbind(data.frame(
    rate = same_benefits_one,
    genotype = '[0, 1]', num_muts = '1',
    scenario = names[1]
  ))%>%
  rbind(data.frame(
    rate = same_benefits_two,
    genotype = '[1, 1]', num_muts = '2',
    scenario = names[1]
  )) %>%
  rbind(data.frame(
    rate = only_one,
    genotype = '[1, 0]', num_muts = '1',
    scenario = names[2]
  )) %>%
  rbind(data.frame(
    rate = low_one,
    genotype = '[0, 1]', num_muts = '1',
    scenario = names[3]
  ))%>%
  rbind(data.frame(
    rate = high_one,
    genotype = '[1, 0]', num_muts = '1',
    scenario = names[3]
  )) %>%
  rbind(data.frame(
    rate = low_high_two,
    genotype = '[1, 1]', num_muts = '2',
    scenario = names[3]
  )) 

nocost_death_additive$concentration = rep(
  concentrations, nrow(nocost_death_additive)/length(concentrations))

nocost_death_additive$cost = 'No cost'

# Create plot
ggplot(nocost_death_additive, aes(x=concentration, y=rate, color = genotype, linetype = num_muts)) +
  facet_grid(scenario~.) +
  geom_line(linewidth = 1.5) +
  theme_bw(base_size = 18) +
  labs(color = 'genotype') +
  ylab('Growth rate (arbitrary units)') +
  # xlim(0.05,100)+
  scale_x_log10(label=label_power10(magnitude_only = TRUE), limits = c(0.5, 1e3)) + 
  # scale_x_log10(limits = c(0.05, 1e3))+
  xlab('Drug concentration (arbitrary units)') 



##### Cost Drug effects on death (not growth) rate, additive (not multiplicative) on MIC #####

cost_vector <- c(0.2, 0.2, 0.2)

# Same benefits, two loci
double_same <- c(high_ben, high_ben, 0)
single_same <- c(high_ben, 0, 0)

# One locus R
single_only <- c(high_ben, 0, 0)

# Low, high benefit, two loci
double_low_high <- c(low_ben, high_ben, 0)
single_low <- c(low_ben, 0, 0)
single_high <- c(0, high_ben, 0)

wildtype <- growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, c(0,0,0), c(0,0,0))

# Calculate growth rates for different scenarios
same_benefits_two <- growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, double_same, c(1,1,0))
same_benefits_one <- growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, single_same, c(1,0,0))

only_one <- growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, single_only, c(1,0,0))

low_high_two <- growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, double_low_high, c(1,1,0))
low_one <- growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, single_low, c(1,0,0))
high_one <- growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, single_high, c(0,1,0))

# Create dataframe (Note: you may need to adjust the variable names)
# The Python code references variables not defined in the notebook

names = c('same/similar benefits \n for beneficial mutations', 
          'only one beneficial \n mutation available',
          'low or high benefit \n for beneficial mutations')

cost_death_additive <- data.frame(
  rate = rep(wildtype, 3),
  genotype = '[0, 0]', num_muts = '0',
  scenario = rep(names, each = length(wildtype))
) %>%
  rbind(data.frame(
    rate = same_benefits_one,
    genotype = '[1, 0]', num_muts = '1',
    scenario = names[1]
  )) %>%
  rbind(data.frame(
    rate = same_benefits_one,
    genotype = '[0, 1]', num_muts = '1',
    scenario = names[1]
  )) %>%
  rbind(data.frame(
    rate = same_benefits_two,
    genotype = '[1, 1]', num_muts = '2',
    scenario = names[1]
  )) %>%
  rbind(data.frame(
    rate = only_one,
    genotype = '[1, 0]', num_muts = '1',
    scenario = names[2]
  )) %>%
  rbind(data.frame(
    rate = low_one,
    genotype = '[0, 1]', num_muts = '1',
    scenario = names[3]
  ))%>%
  rbind(data.frame(
    rate = high_one,
    genotype = '[1, 0]', num_muts = '1',
    scenario = names[3]
  )) %>%
  rbind(data.frame(
    rate = low_high_two,
    genotype = '[1, 1]', num_muts = '2',
    scenario = names[3]
  )) 

cost_death_additive$concentration = rep(
  concentrations, nrow(cost_death_additive)/length(concentrations))

cost_death_additive$cost = 'Cost'

death_additive_cost_nocost = cost_death_additive %>%
  rbind(nocost_death_additive) %>%
  dplyr::mutate(scenario = factor(
    scenario, 
    levels = c('only one beneficial \n mutation available',
               'same/similar benefits \n for beneficial mutations',
               'low or high benefit \n for beneficial mutations')),
    cost = factor(cost, levels = c('No cost', 'Cost')))



# Create plot
ggplot(death_additive_cost_nocost, 
       aes(x=concentration, y=rate, color = genotype, 
           linetype = genotype, size=genotype)) +
  facet_grid(scenario~cost) +
  geom_line() +
  theme_bw(base_size = 23) +
  labs(color = 'genotype')+
  scale_color_manual(values = c('black',hex))+
  scale_size_manual(values = c(0.7, 1.5, 1.5, 2))+
  scale_linetype_manual(values=c("solid", "twodash", "dashed", "dotted")) +
  ylab('Growth rate (arbitrary units)') +
  # xlim(0.05,100)+
  # scale_x_log10(limits = c(0.05, 1e3))+
  xlab('Drug concentration (arbitrary units)') +
  # gtable_add_grob(list(
  #   rectGrob(gp = gpar(col = NA, fill = gray(0.5))),
  #   textGrob("Variable 2", gp = gpar(col = gray(1)))),
  #    3, 4, 3, 6, name = paste(runif(2))) + 
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Genotype - phenotype map", breaks = NULL, labels = NULL)) +
  # scale_x_continuous(sec.axis = sec_axis(~ . , name = "Cost", breaks = NULL, labels = NULL)) 
  scale_x_log10(label=label_power10(magnitude_only = TRUE), limits = c(1, 1e3)) 


###### Repeatable, bidirectional cross R no cost ########

all_genotypes = expand.grid(rep(list(c(0, 1)), 1))

benefits_vector = cbind(c(high_ben),c(high_ben))
cost_vector = c(0)

bidir_rep <- data.frame(
  rate = rep(
    growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, c(0), c(0)), 
    2),
  genotype= '0', num_muts = '0', 
  scenario = rep(c('Drug A', 'Drug B'), each = length(concentrations)) 
)

for (i in 2:nrow(all_genotypes)){
  geno = all_genotypes[i,]
  bidir_rep = rbind(
    bidir_rep, 
    data.frame(
      rate = growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, benefits_vector[,1], geno),
      genotype = paste(geno, collapse = ' '), num_muts = sum(geno!=0), 
      scenario = 'Drug A'
    ),
    data.frame(
      rate = growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, benefits_vector[,2], geno),
      genotype = paste(geno, collapse = ' '), num_muts = sum(geno!=0), 
      scenario = 'Drug B'
    ))
}

bidir_rep$concentration = rep(
  concentrations, nrow(bidir_rep)/length(concentrations))

bidir_rep = bidir_rep %>%
  dplyr::mutate(
    xmin_noMIC_single = 10^0.5,
    xmax_MIC_single = case_when(
      scenario == 'Drug A' ~ 10^(2.34),
      scenario == 'Drug B' ~ 10^(2.34)))

hex <-  palette.colors(palette = "Okabe-Ito")[2]

fix_spacing = 'genotype  '
# Create plot
ggplot(bidir_rep, 
       aes(x=concentration, y=rate, color = genotype, size=genotype)) +
  facet_grid(scenario~.) +
  geom_rect(
    data = unique(bidir_rep[c("xmin_noMIC_single","xmax_MIC_single","scenario")]),
    aes(xmin =xmin_noMIC_single, xmax = xmax_MIC_single,
        ymin=-Inf, ymax=Inf), 
    inherit.aes = FALSE, alpha=0.4, fill='darkgrey') +
  geom_line(alpha=0.7) +
  theme_bw(base_size = 23) +
  labs(color = 'genotype')+
  scale_color_manual(values = c('black',hex))+
  scale_size_manual(values = c(0.7, 1.5))+
  # scale_linetype_manual(values=c("solid", "dashed")) +
  ylab('Growth rate (arbitrary units)') +
  xlab('Drug concentration (arbitrary units)')  +
  scale_x_log10(label=label_power10(magnitude_only = TRUE), limits = c(1, 1e3)) +
  scale_y_continuous(breaks = c(-1,0,1), minor_breaks = c(-0.75, -0.5, -0.25, 0.25, 0.5, 0.75))+
  guides(linetype = guide_legend(override.aes = list(size = 4), keywidth = 4, keyheight = 1),
         colour=guide_legend(keywidth = 4, keyheight = 1))+
  labs(color = fix_spacing, linetype = fix_spacing, size = fix_spacing)

ggsave('/Users/louagema/Documents/Figures/CollaterEffectsConcept/bidirectional_repeatable.png',
       dpi = 300, height = 6, width = 7.5)



###### Non-repeatable, bidirectional cross R no cost ########

print(k)

all_genotypes = expand.grid(rep(list(c(0, 1)), 3))

# -1 and 0.93* are for visualization purposes only
benefits_vector = cbind(c(0.93*high_ben, high_ben, -1),c(-1, high_ben, 0.93*high_ben))
# 0.21 for illustration only
cost_vector = c(0.24, 0.2, 0.24)
# cost_vector = c(0, 0, 0)

bidir_nonrep <- data.frame(
  rate = rep(
    growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, c(0,0,0), c(0,0,0)), 
    2),
  genotype = '0 0 0', num_muts = '0', 
  scenario = rep(c('Drug A', 'Drug B'), each = length(concentrations)) 
)
# all_genotypes = all_genotypes[c(1,2,3,5,4,6,7,8),]
for (i in 2:nrow(all_genotypes)){
  geno = all_genotypes[i,]
  bidir_nonrep = rbind(
    bidir_nonrep, 
    data.frame(
      rate = growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, benefits_vector[,1], geno),
      genotype = paste(geno, collapse = ' '), num_muts = sum(geno!=0), 
      scenario = 'Drug A'
    ),
    data.frame(
      rate = growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, benefits_vector[,2], geno),
      genotype = paste(geno, collapse = ' '), num_muts = sum(geno!=0), 
      scenario = 'Drug B'
    ))
}

bidir_nonrep$concentration = rep(
  concentrations, nrow(bidir_nonrep)/length(concentrations))

bidir_nonrep = bidir_nonrep %>%
  dplyr::mutate(
    xmin_noMIC_single = 10^1,
    xmax_MIC_single = case_when(
      scenario == 'Drug A' ~ 10^(2.09),
      scenario == 'Drug B' ~ 10^(2.09))) %>% 
  dplyr::mutate(genotype = factor(
    genotype, 
    levels = c('0 0 0', '0 0 1', '1 0 0',  '0 1 0',
               '1 1 0', '1 0 1', '0 1 1', '1 1 1')))

hex <- palette.colors(palette = "Okabe-Ito") #hue_pal()(7)

hex[4] = hex[2]
hex[2] = palette.colors(palette = "Okabe-Ito")[4]
# Create plot
ggplot(bidir_nonrep, 
       aes(x=concentration, y=rate, color = genotype, 
           linetype = genotype, size=genotype)) +
  facet_grid(scenario~.) +
  geom_rect(
    data = unique(bidir_nonrep[c("xmin_noMIC_single","xmax_MIC_single","scenario")]),
    aes(xmin =xmin_noMIC_single, xmax = xmax_MIC_single,
        ymin=-Inf, ymax=Inf), 
    inherit.aes = FALSE, alpha=0.4, fill='darkgrey')  +
  geom_line(data=bidir_nonrep %>% dplyr::filter(genotype == '0 0 0'))+
  geom_line(data=bidir_nonrep %>% dplyr::filter(genotype != '0 0 0'), alpha=0.9)+
  theme_bw(base_size = 23)+
  labs(color = 'genotype')+
  scale_color_manual(values = c('black', hex[2:8]))+
  scale_size_manual(values = c(0.5, rep(1.5, 7)))+
  scale_linetype_manual(values=c(rep("solid", 4), rep("dashed", 3), "dotted")) +
  ylab('Growth rate (arbitrary units)') +
  xlab('Drug concentration (arbitrary units)')  +
  scale_x_log10(label=label_power10(magnitude_only = TRUE), limits = c(1, 1e3)) +
  scale_y_continuous(breaks = c(-1,0,1), minor_breaks = c(-0.75, -0.5, -0.25, 0.25, 0.5, 0.75)) +
  guides(linetype = guide_legend(override.aes = list(size = 4), keywidth = 4, keyheight = 1),
         colour=guide_legend(keywidth = 4, keyheight = 1))


ggsave('/Users/louagema/Documents/Figures/CollaterEffectsConcept/bidirectional_nonrepeatable.png',
       dpi = 300, height = 6, width = 7.5)



###### Non-repeatable in one direction repeatable in other, cross R no cost ########



all_genotypes = expand.grid(rep(list(c(0, 1)), 2))

benefits_vector = cbind(c(0.9*high_ben, high_ben),c(-2, high_ben))
cost_vector = c(0.2, 0.2)
cost_vector = c(0, 0)

bidir_nonreprep <- data.frame(
  rate = rep(
    growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, c(0,0), c(0,0)), 
    2),
  genotype = '0 0', num_muts = '0', 
  scenario = rep(c('Drug A', 'Drug B'), each = length(concentrations)) 
)
# all_genotypes = all_genotypes[c(1,2,3,5,4,6,7,8),]
for (i in 2:nrow(all_genotypes)){
  geno = all_genotypes[i,]
  bidir_nonreprep = rbind(
    bidir_nonreprep, 
    data.frame(
      rate = growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, benefits_vector[,1], geno),
      genotype = paste(geno, collapse = ' '), num_muts = sum(geno!=0), 
      scenario = 'Drug A'
    ),
    data.frame(
      rate = growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, benefits_vector[,2], geno),
      genotype = paste(geno, collapse = ' '), num_muts = sum(geno!=0), 
      scenario = 'Drug B'
    ))
}

bidir_nonreprep$concentration = rep(
  concentrations, nrow(bidir_nonreprep)/length(concentrations))

bidir_nonreprep = bidir_nonreprep %>%
  dplyr::mutate(
    xmin_noMIC_single = 10^0.65,
    xmax_MIC_single = case_when(
      scenario == 'Drug A' ~ 10^(1.65),
      scenario == 'Drug B' ~ 10^(2.35)))  %>%
    dplyr::mutate(genotype = factor(
    genotype,
     levels = c('0 0', '0 1', '1 0',  '1 1')))

hex <- palette.colors(palette = "Okabe-Ito") #hue_pal()(7)

# Create plot
ggplot(bidir_nonreprep, 
       aes(x=concentration, y=rate, color = genotype, 
           linetype = genotype, size=genotype)) +
  facet_grid(scenario~.) +
  geom_rect(
    data = unique(bidir_nonreprep[c("xmin_noMIC_single","xmax_MIC_single","scenario")]),
    aes(xmin =xmin_noMIC_single, xmax = xmax_MIC_single,
        ymin=-Inf, ymax=Inf), 
    inherit.aes = FALSE, alpha=0.4, fill='darkgrey')  +
  geom_line(data=bidir_nonreprep %>% dplyr::filter(genotype != '0 0 0'))+
  geom_line(data=bidir_nonreprep %>% dplyr::filter(genotype == '0 0 0'))+
  theme_bw(base_size = 23)+
  labs(color = 'genotype')+
  scale_color_manual(values = c('black', hex[2:4]))+
  scale_size_manual(values = c(1.5, 1.5, 1.5, 1.5))+
  scale_linetype_manual(values=c(rep("solid", 3), rep("dashed", 3), "dotted", 'solid')) +
  ylab('Growth rate (arbitrary units)') +
  xlab('Drug concentration (arbitrary units)')  +
  scale_x_log10(label=label_power10(magnitude_only = TRUE), limits = c(1, 1e3)) +
  scale_y_continuous(breaks = c(-1,0,1), minor_breaks = c(-0.75, -0.5, -0.25, 0.25, 0.5, 0.75)) +
  guides(linetype = guide_legend(override.aes = list(size = 4), keywidth = 4, keyheight = 1),
         colour=guide_legend(keywidth = 4, keyheight = 1))


ggsave('/Users/louagema/Documents/Figures/CollaterEffectsConcept/bidirectional_nonrepeatablerep.png',
       dpi = 300, height = 6, width = 7.5)


###### Unidirectional repeatable cross R cost ########


all_genotypes = expand.grid(rep(list(c(0, 1)), 2))
k=4
benefits_vector = cbind(c(50, 0),c(50, high_ben))

# k=3
# benefits_vector = cbind(c(low_ben, 0),c(low_ben, high_ben))

cost_vector = c(0.2, 0.2)
# cost_vector = c(0, 0, 0)

unidir_rep <- data.frame(
  rate = rep(
    growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, c(0,0), c(0,0)), 
    2),
  genotype = '0 0', num_muts = '0', 
  scenario = rep(c('Drug A', 'Drug B'), each = length(concentrations)) 
)
for (i in 2:nrow(all_genotypes)){
  geno = all_genotypes[i,]
  unidir_rep = rbind(
    unidir_rep, 
    data.frame(
      rate = growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, benefits_vector[,1], geno),
      genotype = paste(geno, collapse = ' '), num_muts = sum(geno!=0), 
      scenario = 'Drug A'
    ),
    data.frame(
      rate = growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, benefits_vector[,2], geno),
      genotype = paste(geno, collapse = ' '), num_muts = sum(geno!=0), 
      scenario = 'Drug B'
    ))
}

unidir_rep$concentration = rep(
  concentrations, nrow(unidir_rep)/length(concentrations))

unidir_rep = unidir_rep %>%
  dplyr::mutate(
    xmin_noMIC_single = case_when(
      scenario == 'Drug A' ~ 10^(1.05),
      scenario == 'Drug B' ~ 10^(1.85)),
    xmax_MIC_single = case_when(
      scenario == 'Drug A' ~ 10^(1.85),
      scenario == 'Drug B' ~ 10^(2.2)),
    vlines = case_when(
      scenario == 'Drug A' ~ NA_real_,
      scenario == 'Drug B' ~ 10^(1.45)
    )
    )

# unidir_rep = unidir_rep %>%
#   dplyr::mutate(
#     xmin_noMIC_single = case_when(
#       scenario == 'Drug A' ~ 10^(1),
#       scenario == 'Drug B' ~ 10^(2)),
#     xmax_MIC_single = case_when(
#       scenario == 'Drug A' ~ 10^(2),
#       scenario == 'Drug B' ~ 10^(2.2)),
#     vlines = case_when(
#       scenario == 'Drug A' ~ NA_real_,
#       scenario == 'Drug B' ~ 10^(1.45)
#     )
#   )

hex <- palette.colors(palette = "Okabe-Ito") #hue_pal()(7)

fix_spacing = 'genotype  '
# Create plot
ggplot(unidir_rep, 
       aes(x=concentration, y=rate, color = genotype, 
           linetype = genotype, size=genotype)) +
  facet_grid(scenario~.) +
  geom_hline(yintercept = 0, color = 'darkgrey', size=1)+
  geom_rect(
    data = unique(unidir_rep[c("xmin_noMIC_single","xmax_MIC_single","scenario")]),
    aes(xmin =xmin_noMIC_single, xmax = xmax_MIC_single,
        ymin=-Inf, ymax=Inf), 
    inherit.aes = FALSE, alpha=0.4, fill='darkgrey')+
  geom_rect(
    data = unique(unidir_rep[c("xmin_noMIC_single","vlines","scenario")]),
    aes(xmin =vlines, xmax = xmin_noMIC_single,
        ymin=-Inf, ymax=Inf), 
    inherit.aes = FALSE, alpha=0.4, fill='orange') +
  geom_line(data=unidir_rep %>% dplyr::filter(genotype != '0 0'))+
  geom_line(data=unidir_rep %>% dplyr::filter(genotype == '0 0')) +
  theme_bw(base_size = 23)+
  labs(color = 'genotype')+
  scale_color_manual(values = c('black',hex[2:10]))+
  scale_size_manual(values = c(0.5, 1.7, rep(1.5, 7)))+
  scale_linetype_manual(values=c("solid", rep("solid", 2), rep("dashed", 1), "dotted")) +
  ylab('Growth rate (arbitrary units)') +
  xlab('Drug concentration (arbitrary units)')  +
  scale_x_log10(label=label_power10(magnitude_only = TRUE), limits = c(3.3, 1e3)) +
  scale_y_continuous(breaks = c(-1,0,1), minor_breaks = c(-0.75, -0.5, -0.25, 0.25, 0.5, 0.75)) +
  guides(linetype = guide_legend(override.aes = list(size = 4), keywidth = 4, keyheight = 1),
         colour=guide_legend(keywidth = 4, keyheight = 1)) +
  labs(color = fix_spacing, linetype = fix_spacing, size = fix_spacing)


ggsave('/Users/louagema/Documents/Figures/CollaterEffectsConcept/unidirectional_repeatable.png',
       dpi = 300, height = 6, width = 7.5)



###### Unidirectional repeatable cross R cost less extinction ########

MIC_uni = 200
concentrations_uni = 10^seq(1, 4, length.out = 100)
all_genotypes = expand.grid(rep(list(c(0, 1)), 2))
k_uni=3
benefits_vector = cbind(c(50, 0),c(50, 90))
cost_vector = c(0.05, 0.05)
# cost_vector = c(0, 0, 0)

unidir_rep <- data.frame(
  rate = rep(
    growth_rate_genotype(phi_max, phi_min, concentrations_uni, MIC_uni, k_uni, cost_vector, c(0,0), c(0,0)), 
    2),
  genotype = '0 0', num_muts = '0', 
  scenario = rep(c('Drug A', 'Drug B'), each = length(concentrations_uni)) 
)
for (i in 2:nrow(all_genotypes)){
  geno = all_genotypes[i,]
  unidir_rep = rbind(
    unidir_rep, 
    data.frame(
      rate = growth_rate_genotype(phi_max, phi_min, concentrations_uni, MIC_uni, k_uni, cost_vector, benefits_vector[,1], geno),
      genotype = paste(geno, collapse = ' '), num_muts = sum(geno!=0), 
      scenario = 'Drug A'
    ),
    data.frame(
      rate = growth_rate_genotype(phi_max, phi_min, concentrations_uni, MIC_uni, k_uni, cost_vector, benefits_vector[,2], geno),
      genotype = paste(geno, collapse = ' '), num_muts = sum(geno!=0), 
      scenario = 'Drug B'
    ))
}

unidir_rep$concentration = rep(
  concentrations_uni, nrow(unidir_rep)/length(concentrations_uni))

unidir_rep = unidir_rep %>%
  dplyr::mutate(
    xmin_noMIC_single = case_when(
      scenario == 'Drug A' ~ 10^(1.05),
      scenario == 'Drug B' ~ 10^(1.85)),
    xmax_MIC_single = case_when(
      scenario == 'Drug A' ~ 10^(1.85),
      scenario == 'Drug B' ~ 10^(2.2)),
    vlines = case_when(
      scenario == 'Drug A' ~ NA_real_,
      scenario == 'Drug B' ~ 10^(1.45)
    )
  )

hex <- palette.colors(palette = "Okabe-Ito") #hue_pal()(7)

fix_spacing = 'genotype  '
# Create plot
ggplot(unidir_rep, 
       aes(x=concentration, y=rate, color = genotype, 
           linetype = genotype, size=genotype)) +
  facet_grid(scenario~.) +
  geom_rect(
    data = unique(unidir_rep[c("xmin_noMIC_single","xmax_MIC_single","scenario")]),
    aes(xmin =xmin_noMIC_single, xmax = xmax_MIC_single,
        ymin=-Inf, ymax=Inf), 
    inherit.aes = FALSE, alpha=0.4, fill='darkgrey')+
  geom_rect(
    data = unique(unidir_rep[c("xmin_noMIC_single","vlines","scenario")]),
    aes(xmin =vlines, xmax = xmin_noMIC_single,
        ymin=-Inf, ymax=Inf), 
    inherit.aes = FALSE, alpha=0.4, fill='orange') +
  geom_line(data=unidir_rep %>% dplyr::filter(genotype != '0 0'))+
  geom_line(data=unidir_rep %>% dplyr::filter(genotype == '0 0')) +
  theme_bw(base_size = 23)+
  labs(color = 'genotype')+
  scale_color_manual(values = c('black',hex[2:10]))+
  scale_size_manual(values = c(0.5, 1.7, rep(1.5, 7)))+
  scale_linetype_manual(values=c("solid", rep("solid", 2), rep("dashed", 1), "dotted")) +
  ylab('Growth rate (arbitrary units)') +
  xlab('Drug concentration (arbitrary units)')  +
  scale_x_log10(label=label_power10(magnitude_only = TRUE), limits = c(3.3, 1e3)) +
  scale_y_continuous(breaks = c(-1,0,1), minor_breaks = c(-0.75, -0.5, -0.25, 0.25, 0.5, 0.75)) +
  guides(linetype = guide_legend(override.aes = list(size = 4), keywidth = 4, keyheight = 1),
         colour=guide_legend(keywidth = 4, keyheight = 1)) +
  labs(color = fix_spacing, linetype = fix_spacing, size = fix_spacing)


ggsave('/Users/louagema/Documents/Figures/CollaterEffectsConcept/unidirectional_repeatable.png',
       dpi = 300, height = 6, width = 7.5)



###### Collateral sensitivity of cross resistance non-repeatable ######

benefits_vector = cbind(c(high_ben,high_ben,low_ben),c(0,-MIC/1.5,low_ben))
cost_vector <- c(0.2, 0.2, 0.2)

all_genotypes = expand.grid(rep(list(c(0, 1)), 3))

CSCR_nocost_death_additive <- data.frame(
  rate = rep(
    growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, c(0,0,0), c(0,0,0)), 
    2),
  genotype = '0 0 0', num_muts = '0', 
  scenario = rep(c('Drug A', 'Drug B'), each = length(concentrations)) 
)
for (i in 2:nrow(all_genotypes)){
  geno = all_genotypes[i,]
  CSCR_nocost_death_additive = rbind(
    CSCR_nocost_death_additive, 
    data.frame(
      rate = growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, benefits_vector[,1], geno),
      genotype = paste(geno, collapse = ' '), num_muts = sum(geno!=0), 
      scenario = 'Drug A'
    ),
    data.frame(
      rate = growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, benefits_vector[,2], geno),
      genotype = paste(geno, collapse = ' '), num_muts = sum(geno!=0), 
      scenario = 'Drug B'
    ))
}

CSCR_nocost_death_additive$concentration = rep(
  concentrations, nrow(CSCR_nocost_death_additive)/length(concentrations))

CSCR_nocost_death_additive = CSCR_nocost_death_additive %>%
  dplyr::mutate(
    xmin_noMIC_single = 10^0.5,
    xmax_noMIC_single = case_when(
      scenario == 'Drug A' ~ 10,
      scenario == 'Drug B' ~ Inf),
    xmax_MIC_single = case_when(
      scenario == 'Drug A' ~ 10,
      scenario == 'Drug B' ~ 10^1.7))


hex <- hue_pal()(8)

# Create plot
ggplot(CSCR_nocost_death_additive, 
       aes(x=concentration, y=rate, color = genotype, linetype = genotype)) +
  facet_grid(~scenario) +
  geom_line(linewidth=1.5) +
  theme_bw(base_size = 23) +
  labs(color = 'genotype')+
  scale_color_manual(values = c('black',hex))+
  # scale_size_manual(values = c(0.7, 1.5, 1.5, 2))+
  scale_linetype_manual(values=c(
    "solid", "twodash", "twodash", "twodash", 
    "dashed", "dashed", "dashed", "dotted")) +
  ylab('Growth rate (arbitrary units)') +
  xlab('Drug concentration (arbitrary units)') +
  geom_vline(xintercept = c(10,100), alpha=0.6)+
  scale_x_log10(label=label_power10(magnitude_only = TRUE), limits = c(1, 1e3)) 


# Create plot
ggplot(CSCR_nocost_death_additive, 
       aes(x=concentration, y=rate, color = genotype, 
           linetype = genotype, size=genotype)) +
  facet_grid(scenario~.) +
  annotate('rect', xmin=10^0.5, xmax=10, ymin=-Inf, ymax=Inf, alpha=.7, fill='lightgrey') +
  # annotate('rect', xmin=10, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=.2, fill='pink') +
  geom_line() +
  theme_bw(base_size = 23) +
  labs(color = 'genotype')+
  scale_color_manual(values = c('black',hex))+
  scale_size_manual(values = c(0.7, 1.5, 1.5, 2.3))+
  scale_linetype_manual(values=c("solid", "twodash", "dashed", "dotted")) +
  ylab('Growth rate (arbitrary units)') +
  xlab('Drug concentration (arbitrary units)')  +
  scale_x_log10(label=label_power10(magnitude_only = TRUE), limits = c(1, 1e3)) 

###### Collateral sensitivity of cross resistance repeatable ######

k=3
benefits_vector = cbind(c(high_ben,low_ben),c(-MIC/1.5,low_ben))
cost_vector <- c(0.2, 0.2)

all_genotypes = expand.grid(rep(list(c(0, 1)), 2))

CSCR_nocost_death_additive <- data.frame(
  rate = rep(
    growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, c(0,0), c(0,0)), 
    2),
  genotype = '0 0', num_muts = 0, 
  scenario = rep(c('Drug A', 'Drug B'), each = length(concentrations)) 
)
for (i in 2:nrow(all_genotypes)){
  geno = all_genotypes[i,]
  CSCR_nocost_death_additive = rbind(
    CSCR_nocost_death_additive, 
    data.frame(
      rate = growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, benefits_vector[,1], geno),
      genotype = paste(geno, collapse = ' '), num_muts = sum(geno!=0), 
      scenario = 'Drug A'
    ),
    data.frame(
      rate = growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, benefits_vector[,2], geno),
      genotype = paste(geno, collapse = ' '), num_muts = sum(geno!=0), 
      scenario = 'Drug B'
    ))
}

CSCR_nocost_death_additive$concentration = rep(
  concentrations, nrow(CSCR_nocost_death_additive)/length(concentrations))

CSCR_nocost_death_additive = CSCR_nocost_death_additive %>%
  dplyr::mutate(
    xmin_single = case_when(
      scenario == 'Drug A' ~ 10^(1+0.35),
      scenario == 'Drug B' ~ 10^(1)),
    xmax_MIC_single = case_when(
      scenario == 'Drug A' ~ 10^(2.1),
      scenario == 'Drug B' ~ 10^(2)),
    vlines = case_when(
      scenario == 'Drug A' ~  10^(1),
      scenario == 'Drug B' ~ NA_real_
    ))

hex <- palette.colors(palette = "Okabe-Ito") #hue_pal()(7)

fix_spacing = 'genotype  '
# Create plot
ggplot(CSCR_nocost_death_additive, 
       aes(x=concentration, y=rate, color = genotype, 
           linetype = genotype, size=genotype)) +
  facet_grid(scenario~.) +
  geom_rect(
    data = unique(CSCR_nocost_death_additive[c("xmin_single","xmax_MIC_single","scenario")]),
    aes(xmin =xmin_single, xmax = xmax_MIC_single,
                ymin=-Inf, ymax=Inf), 
            inherit.aes = FALSE, alpha=0.4, fill='darkgrey') +
  geom_rect(
    data = unique(CSCR_nocost_death_additive[c("xmin_single","vlines","scenario")]),
    aes(xmin =vlines, xmax = xmin_single,
        ymin=-Inf, ymax=Inf), 
    inherit.aes = FALSE, alpha=0.4, fill='orange') +
  # geom_vline(aes(xintercept = vlines), color = 'darkgrey', 
  #            linetype = 'dashed', linewidth = 0.8)+
  geom_line() +
  theme_bw(base_size = 23) +
  labs(color = 'genotype')+
  scale_color_manual(values = c('black',hex[2:5]))+
  scale_size_manual(values = c(0.7, 1.5, 1.5, 1.5))+
  scale_linetype_manual(values=c("solid", "solid", "solid", "dashed")) +
  ylab('Growth rate (arbitrary units)') +
  xlab('Drug concentration (arbitrary units)')  +
  scale_x_log10(label=label_power10(magnitude_only = TRUE), limits = c(1, 1e3))+
  scale_y_continuous(breaks = c(-1,0,1), minor_breaks = c(-0.75, -0.5, -0.25, 0.25, 0.5, 0.75)) +
  guides(linetype = guide_legend(override.aes = list(size = 4), keywidth = 4, keyheight = 1),
         colour=guide_legend(keywidth = 4, keyheight = 1))+
  labs(color = fix_spacing, linetype = fix_spacing, size = fix_spacing)

ggsave('/Users/louagema/Documents/Figures/CollaterEffectsConcept/CRCS_repeatable.png',
       dpi = 300, height = 6, width = 7.5)



###### Collateral sensitivity of cross resistance repeatable ######



k=3
benefits_vector = cbind(c(high_ben, -MIC/1.5),c(-MIC/1.5,high_ben))
cost_vector <- c(0.2, 0.2)

all_genotypes = expand.grid(rep(list(c(0, 1)), 2))

CSCR_nocost_death_additive <- data.frame(
  rate = rep(
    growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, c(0,0), c(0,0)), 
    2),
  genotype = '0 0', num_muts = 0, 
  scenario = rep(c('Drug A', 'Drug B'), each = length(concentrations)) 
)
for (i in 2:nrow(all_genotypes)){
  geno = all_genotypes[i,]
  CSCR_nocost_death_additive = rbind(
    CSCR_nocost_death_additive, 
    data.frame(
      rate = growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, benefits_vector[,1], geno),
      genotype = paste(geno, collapse = ' '), num_muts = sum(geno!=0), 
      scenario = 'Drug A'
    ),
    data.frame(
      rate = growth_rate_genotype(phi_max, phi_min, concentrations, MIC, k, cost_vector, benefits_vector[,2], geno),
      genotype = paste(geno, collapse = ' '), num_muts = sum(geno!=0), 
      scenario = 'Drug B'
    ))
}

CSCR_nocost_death_additive$concentration = rep(
  concentrations, nrow(CSCR_nocost_death_additive)/length(concentrations))

CSCR_nocost_death_additive = CSCR_nocost_death_additive %>%
  dplyr::mutate(
    xmin_single = case_when(
      scenario == 'Drug A' ~ 10^(1),
      scenario == 'Drug B' ~ 10^(1)),
    xmax_MIC_single = case_when(
      scenario == 'Drug A' ~ 10^(2.35),
      scenario == 'Drug B' ~ 10^(2.35)),
    vlines = case_when(
      scenario == 'Drug A' ~  NA_real_,
      scenario == 'Drug B' ~ NA_real_
    ))

hex <- palette.colors(palette = "Okabe-Ito") #hue_pal()(7)

fix_spacing = 'genotype  '
# Create plot
ggplot(CSCR_nocost_death_additive, 
       aes(x=concentration, y=rate, color = genotype, 
           linetype = genotype, size=genotype)) +
  facet_grid(scenario~.) +
  geom_rect(
    data = unique(CSCR_nocost_death_additive[c("xmin_single","xmax_MIC_single","scenario")]),
    aes(xmin =xmin_single, xmax = xmax_MIC_single,
        ymin=-Inf, ymax=Inf), 
    inherit.aes = FALSE, alpha=0.4, fill='darkgrey') +
  geom_rect(
    data = unique(CSCR_nocost_death_additive[c("xmin_single","vlines","scenario")]),
    aes(xmin =vlines, xmax = xmin_single,
        ymin=-Inf, ymax=Inf), 
    inherit.aes = FALSE, alpha=0.4, fill='orange') +
  # geom_vline(aes(xintercept = vlines), color = 'darkgrey', 
  #            linetype = 'dashed', linewidth = 0.8)+
  geom_line() +
  theme_bw(base_size = 23) +
  labs(color = 'genotype')+
  scale_color_manual(values = c('black',hex[2:5]))+
  scale_size_manual(values = c(0.7, 1.5, 1.5, 1.5))+
  scale_linetype_manual(values=c("solid", "solid", "solid", "dashed")) +
  ylab('Growth rate (arbitrary units)') +
  xlab('Drug concentration (arbitrary units)')  +
  scale_x_log10(label=label_power10(magnitude_only = TRUE), limits = c(1, 1e3))+
  scale_y_continuous(breaks = c(-1,0,1), minor_breaks = c(-0.75, -0.5, -0.25, 0.25, 0.5, 0.75)) +
  guides(linetype = guide_legend(override.aes = list(size = 4), keywidth = 4, keyheight = 1),
         colour=guide_legend(keywidth = 4, keyheight = 1))+
  labs(color = fix_spacing, linetype = fix_spacing, size = fix_spacing)

ggsave('/Users/louagema/Documents/Figures/CollaterEffectsConcept/CS_bidir_repeatable.png',
       dpi = 300, height = 6, width = 7.5)




#### Biology 26 #####


# wildtype drug A
ggplot(unidir_rep %>% 
         dplyr::filter(scenario =='Drug A') %>%
         dplyr::filter(genotype=='0 0'), 
       aes(x=concentration, y=rate, color = genotype, 
           linetype = genotype, size=genotype)) +
  facet_grid(scenario~.) +
  geom_hline(yintercept = 0, color = 'darkgrey', size=1)+
  geom_line()+
  theme_bw(base_size = 23)+
  labs(color = 'genotype')+
  scale_color_manual(values = c('black',hex[3:10]))+
  scale_size_manual(values = c(0.5, 1.7, rep(1.5, 7)))+
  scale_linetype_manual(values=c("solid", rep("solid", 2), rep("dashed", 1), "dotted")) +
  ylab('Growth rate (AU)') +
  xlab('Drug concentration (AU)')  +
  scale_x_log10(label=label_power10(magnitude_only = TRUE), limits = c(3.3, 1e3)) +
  scale_y_continuous(breaks = c(-1,0,1), minor_breaks = c(-0.75, -0.5, -0.25, 0.25, 0.5, 0.75)) +
  guides(linetype = guide_legend(override.aes = list(size = 4), keywidth = 4, keyheight = 1),
         colour=guide_legend(keywidth = 4, keyheight = 1)) +
  labs(color = fix_spacing, linetype = fix_spacing, size = fix_spacing)

# wildtype and single mutant drug A


ggplot(unidir_rep %>% 
         dplyr::filter(scenario =='Drug A') %>%
         dplyr::filter(genotype %in% c('0 0', '1 0')), 
       aes(x=concentration, y=rate, color = genotype, 
           linetype = genotype, size=genotype)) +
  facet_grid(scenario~.) +
  geom_hline(yintercept = 0, color = 'darkgrey', size=1)+
  geom_rect(
    data = unique(unidir_rep %>% 
                    dplyr::filter(scenario =='Drug A') %>%
                    dplyr::select(c("xmin_noMIC_single","xmax_MIC_single","scenario"))),
    aes(xmin =xmin_noMIC_single, xmax = xmax_MIC_single,
        ymin=-Inf, ymax=Inf), 
    inherit.aes = FALSE, alpha=0.4, fill='darkgrey')+
  geom_line()+
  theme_bw(base_size = 23)+
  labs(color = 'genotype')+
  scale_color_manual(values = c('black',hex[3:10]))+
  scale_size_manual(values = c(0.5, 1.7, rep(1.5, 7)))+
  scale_linetype_manual(values=c("solid", rep("solid", 2), rep("dashed", 1), "dotted")) +
  ylab('Growth rate (AU)') +
  xlab('Drug concentration (AU)')  +
  scale_x_log10(label=label_power10(magnitude_only = TRUE), limits = c(3.3, 1e3)) +
  scale_y_continuous(breaks = c(-1,0,1), minor_breaks = c(-0.75, -0.5, -0.25, 0.25, 0.5, 0.75)) +
  guides(linetype = guide_legend(override.aes = list(size = 4), keywidth = 4, keyheight = 1),
         colour=guide_legend(keywidth = 4, keyheight = 1)) +
  labs(color = fix_spacing, linetype = fix_spacing, size = fix_spacing)

# wildtype and single mutant drug A and B

unidir_rep_mod = unidir_rep %>%
  dplyr::mutate(xmin_noMIC_single = 
                  case_when(scenario=='Drug B' ~ 11.22,
                            .default = xmin_noMIC_single),
                xmax_MIC_single = 
                  case_when(scenario=='Drug B' ~ 70.79,
                            .default = xmax_MIC_single),
                xmax_orange = 
                  case_when(scenario=='Drug B' ~ 70.79,
                            .default = NA_real_))
ggplot(unidir_rep_mod %>%
         dplyr::filter(genotype %in% c('0 0', '1 0')), 
       aes(x=concentration, y=rate, color = genotype, 
           linetype = genotype, size=genotype)) +
  facet_grid(scenario~.) +
  geom_hline(yintercept = 0, color = 'darkgrey', size=1)+
  geom_rect(
    data = unique(unidir_rep_mod %>%
                    dplyr::select(c("xmin_noMIC_single","xmax_MIC_single","scenario"))),
    aes(xmin =xmin_noMIC_single, xmax = xmax_MIC_single,
        ymin=-Inf, ymax=Inf), 
    inherit.aes = FALSE, alpha=0.4, fill='darkgrey')+
  geom_line()+
  theme_bw(base_size = 23)+
  labs(color = 'genotype')+
  scale_color_manual(values = c('black',hex[3:10]))+
  scale_size_manual(values = c(0.5, 1.7, rep(1.5, 7)))+
  scale_linetype_manual(values=c("solid", rep("solid", 2), rep("dashed", 1), "dotted")) +
  ylab('Growth rate (AU)') +
  xlab('Drug concentration (AU)')  +
  scale_x_log10(label=label_power10(magnitude_only = TRUE), limits = c(3.3, 1e3)) +
  scale_y_continuous(breaks = c(-1,0,1), minor_breaks = c(-0.75, -0.5, -0.25, 0.25, 0.5, 0.75)) +
  guides(linetype = guide_legend(override.aes = list(size = 4), keywidth = 4, keyheight = 1),
         colour=guide_legend(keywidth = 4, keyheight = 1)) +
  labs(color = fix_spacing, linetype = fix_spacing, size = fix_spacing)

ggsave("Biology26_bidirectional_repeatable.png", width = 7.65, height = 5.7)


# non-reciprocity
ggplot(unidir_rep, 
       aes(x=concentration, y=rate, color = genotype, 
           linetype = genotype, size=genotype)) +
  facet_grid(scenario~.) +
  geom_hline(yintercept = 0, color = 'darkgrey', size=1)+
  geom_rect(
    data = unique(unidir_rep %>% 
                    dplyr::filter(grepl('Drug A', scenario)) %>%
                    dplyr::select(c("xmin_noMIC_single","xmax_MIC_single","scenario"))),
    aes(xmin =xmin_noMIC_single, xmax = xmax_MIC_single,
        ymin=-Inf, ymax=Inf), 
    inherit.aes = FALSE, alpha=0.4, fill='darkgrey')+
  geom_rect(
    data = unique(unidir_rep %>% 
                    dplyr::filter(grepl('Drug B', scenario)) %>%
                    dplyr::select(c("xmin_noMIC_single","xmax_MIC_single","scenario"))),
    aes(xmin =xmin_noMIC_single, xmax = xmax_MIC_single,
        ymin=-Inf, ymax=Inf), 
    inherit.aes = FALSE, alpha=0.3, fill='blue')+
  geom_rect(
    data = unique(unidir_rep[c("xmin_noMIC_single","vlines","scenario")]),
    aes(xmin =vlines, xmax = xmin_noMIC_single,
        ymin=-Inf, ymax=Inf), 
    inherit.aes = FALSE, alpha=0.3, fill='red') +
  geom_line()+
  theme_bw(base_size = 23)+
  labs(color = 'genotype')+
  scale_color_manual(values = c('black',hex[2:10]))+
  scale_size_manual(values = c(0.5, 1.7, rep(1.5, 7)))+
  scale_linetype_manual(values=c("solid", rep("solid", 2), rep("dashed", 1), "dotted")) +
  ylab('Growth rate (AU)') +
  xlab('Drug concentration (AU)')  +
  scale_x_log10(label=label_power10(magnitude_only = TRUE), limits = c(3.3, 1e3)) +
  scale_y_continuous(breaks = c(-1,0,1), minor_breaks = c(-0.75, -0.5, -0.25, 0.25, 0.5, 0.75)) +
  guides(linetype = guide_legend(override.aes = list(size = 4), keywidth = 4, keyheight = 1),
         colour=guide_legend(keywidth = 4, keyheight = 1)) +
  labs(color = fix_spacing, linetype = fix_spacing, size = fix_spacing)

ggsave("Biology26_unidirectional_repeatable.png", width = 7.65, height = 5.7)

# non-repeatability
unidir_rep_mod = unidir_rep %>%
  dplyr::mutate(xmin_noMIC_single = 
                  case_when(scenario=='Drug B' ~ 11.22,
                            .default = xmin_noMIC_single),
                xmax_MIC_single = 
                  case_when(scenario=='Drug B' ~ 28.18,
                            .default = xmax_MIC_single),
                xmax_orange = 
                  case_when(scenario=='Drug B' ~ 70.79,
                            .default = NA_real_))

ggplot(unidir_rep_mod, 
       aes(x=concentration, y=rate, color = genotype, 
           linetype = genotype, size=genotype)) +
  facet_grid(scenario~.) +
  geom_hline(yintercept = 0, color = 'darkgrey', size=1)+
  geom_rect(
    data = unique(unidir_rep_mod %>% 
                    dplyr::filter(grepl('Drug A', scenario)) %>%
                    dplyr::select(c("xmin_noMIC_single","xmax_MIC_single","scenario"))),
    aes(xmin =xmin_noMIC_single, xmax = xmax_MIC_single,
        ymin=-Inf, ymax=Inf), 
    inherit.aes = FALSE, alpha=0.4, fill='darkgrey')+
  geom_rect(
    data = unique(unidir_rep_mod %>% 
                    dplyr::filter(grepl('Drug B', scenario)) %>%
                    dplyr::select(c("xmin_noMIC_single","xmax_MIC_single","scenario"))),
    aes(xmin =xmin_noMIC_single, xmax = xmax_MIC_single,
        ymin=-Inf, ymax=Inf), 
    inherit.aes = FALSE, alpha=0.3, fill='green')+
  geom_rect(
    data = unique(unidir_rep_mod[c("xmin_noMIC_single","vlines","scenario", 'xmax_orange')]),
    aes(xmin =vlines, xmax = xmax_orange,
        ymin=-Inf, ymax=Inf), 
    inherit.aes = FALSE, alpha=0.3, fill='green') +
  geom_line()+
  theme_bw(base_size = 23)+
  labs(color = 'genotype')+
  scale_color_manual(values = c('black',hex[2:10]))+
  scale_size_manual(values = c(0.5, 1.7, rep(1.5, 7)))+
  scale_linetype_manual(values=c("solid", rep("solid", 2), rep("dashed", 1), "dotted")) +
  ylab('Growth rate (AU)') +
  xlab('Drug concentration (AU)')  +
  scale_x_log10(label=label_power10(magnitude_only = TRUE), limits = c(3.3, 1e3)) +
  scale_y_continuous(breaks = c(-1,0,1), minor_breaks = c(-0.75, -0.5, -0.25, 0.25, 0.5, 0.75)) +
  guides(linetype = guide_legend(override.aes = list(size = 4), keywidth = 4, keyheight = 1),
         colour=guide_legend(keywidth = 4, keyheight = 1)) +
  labs(color = fix_spacing, linetype = fix_spacing, size = fix_spacing)

ggsave("Biology26_bidirectional_nonrepeatable.png", width = 7.65, height = 5.7)


