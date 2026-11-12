library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(do)
library(ggbeeswarm)
library(ggridges)
library(ggpubr)
library(stringr)
library(ggunchull)
library(purrr)
library(ggrepel)
library(ggbreak)
library(latex2exp)

#### Comments #####
# !!! currently drug names do not match between papers
# cannot compare grouped data across papers
# within paper definitiono f repeatability, directionality is fine
# and can look at that as long as group per paper per drug pair. If exclude per paper, is wrong
# e.g. CEF in Barbosa and Sakenova
# STRUCTURE

# most data processing is done in functions except for initial formatting and adding study details


# To include datasets
#   * Liu 2023
#   * Lázár 2014 
#   * Lazar 2013
#   * > Jahn 2017? Different experimental methods, collateral effects
#   * Kavanaugh 2020 / Flanagan 2020 same dataset?
#   * Rodriguez de Evgrafov 2020
#   * Barbosa 2019?
#   * Wang 2025
#   * Li 2026
#   * Ma 2023
#   * Suzuki 2014
#   * colclough 2019 seaweed and antibiotics mgm

# Not sure if they have data
#   * Chevereau 2015 also knockouts?
#   * Barbosa 2017
#   * Zhao 2016
#   * All searches of Zotero "collateral sensitivity" of 10/05/2025 and earlier
#   18. Lázár, V. et al. Genome-wide analysis captures the determinants
#       of the antibiotic cross-resistance interaction network.
#      Nat. Commun. 5, 4352 (2014).
#   19. Oz, T. et al. Strength of selection pressure is an important
#       parameter contributing to the complexity of antibiotic resistance
#       evolution. Mol. Biol. Evol. 31, 2387–2401 (2014).
#   * Aulin 2021 references, see highlights in goodnotes
#   * Lazar 2013, per repeat CS binarization? And any MIC determination or only relative fitness and at what drug concentration???

# Not useful data
#   * Chauhan 2025 constructed strains themselves and only evolved in CIP at various conc gradients
#   * review? Pal et al 2015 similar to Imamovic 2013
#   * clinical ATLAS database Herencias 2024, unclear how can integrate as can mean drugs are commonly administered together and not collateral effect of evolution to one drug
#           can check for cross R of coll S profiles in in vitro data
#           then see how often these pairs have the same profile BUT if not reciprocal but unidirectional then depends on drug exposure (drug A or B)?. 
#       if ATLAS then also CANWARD
#   * Lazar 2013 no MIC, does have relative fitness to compare cost to benefit?
#         show chemical similarity same as Barbosa
#   * Ardell, not primary data
#   * Santos Lopez 2019: 
#         only evolved in one drug, cannot find raw data. 
#         Biofilm vs planktonic collateral profiles 
#         on different days post selection start
#   * James 2014
#         Only one direction
#   * Nichol 2019
#         Only one direction
#   * Munck Sommer 
#         Combinations or cycles
#   * Oz 2014
#   * CHeck Roemhild
#   * Horinouchi non antibiotic stressors
#   * Liu 2023 not sure what evolved against in data


# Not available data
#   * CANWARD study Karlowski 2019 E. coli (CANNOT FIND IT :((( )
#   * Barbosa 2017: nice heatmap fig 3 pseudomonas
#   * Imamovic 2013, had to manually extract, also did they round?

# Not useful yet
#   * Chevereau, deletion library. Could compare effect size and direction but is seperate thing
#   * Dunne 2022 enterobacterales. Only show Susceptible or R for different species and drugs so no coll S can be derived only co-occurrence of R due to whatever reason. Not MIC. this I think https://bmcinfectdis.biomedcentral.com/articles/10.1186/s12879-022-07167-y
#   * Maltas 2025
#         many days post selection, evaluate collateral profile, I did not include all days. Done for 5 drugs
#           frequency of effect depending on time of selection (fig3)
#   * Barbosa 2017
#         * quantifies cost yield, growth rate, lag time
#   * Nichol 2019
#         Only one direction but 60 replicates, cefotaxime evolution tested against 9 drugs
#   * Chowdhury 2025: barely any CS, fig. 3
#         Experiment only one direction
#         Clinical dataset, use correlation as proxy for cross R, coll S.
#   * Cherny 
#   * Liu 2023, single mutant library (from evolved in one drug or MDR plasmids) have data

# Open questions about data
# Imamovic:
#   * Repeats?: if did, how many had to show effect to show, showed mean or median or?
#   * Rounding?: x2 = 2<x<4?
#   * Show different strains in suppl??
#   * Clinical: seems to show stronger collateral sensitivity
#   * strains evolved resistance to 21 drugs but I have 23 drugs???
# Waller:
#   * Not showing repeats in heatmap means did not evolve resistance to drug evolved in?
# Sakenova: 
#   * overlap between neutral and CS MIC fold change ??? How do they defined Validated interaction?
#   * never tests MIC :( against drug evolved in maybe in other suppl of same paper
#   * !!! add n_repeats!! SHOW SOMETHING, not all repeats though as supposedly 12 lineages Repeats of seperate lineages not shown :(
#   * Evolved in two drugs data
#   * unclear how frequency is defined
#   * Why don't use all drug x drug combinations for evolved x tested?
# Maltas:
#   * Sometimes_collaterally_S_R against CRO why so often?
#   * Why such low MIC fold changes???: NIT has low MIC in heatmap, not in data, CIP has high in heatmap, low in data
# Barbosa 2017:
#   * CON drug not in paper, why? What is it?

# Notes on papers
#   * Maltas 2025 In previous work, we chose this threshold (3-fold for R or S) because it minimized the number of false positives (e.g., all measurements of IC50 in the ancestor strains, across all drugs, fell within this ± 3σWT window) [41,45].

# Summary metrics
#   * How to normalize by number of strains evolved 
#   * Or show per unique drug pair what result is

# TO DO
#   * Show one direction vs other direction effect size = MIC fold change 
#   * Maltas 2018 compare their clustering by drug type and collateral profile with other datasets
#         * Check if within types where expect sensitivity or no effect, cross R is possible?
#         * Imamovic 2013 "These findings indicate that drugs belonging
# to a specific class do not always
# induce the same collateral sensitivity. Instead,
# drug-specific resistance evolution
# plays a role in determining the collateral
# sensitivity profiles."
#   * ATLAS database Herencias 2024, unclear how can integrate as 
#         can mean drugs are commonly administered together and not collateral effect of evolution to one drug
#         can check for cross R of coll S profiles in in vitro data
#         then see how often these pairs have the same profile BUT if not reciprocal but unidirectional then depends on drug exposure (drug A or B)?. 
#         and how define collateral S?
#   * Column for method of evolution
#         * time between drug increases, 
#         * fold change of drug increases
#         * max drug concentration (in terms of MIC fold change relative to wildtype)
#         * soft agar gradient evolution or ...
#         * resolution of MIC fold change relative to wildtype (e.g. MIC determined by two fold dilutions)
#   * Venn diagram UpSetR from drugs tested/evolved in all datasets
#   * Per drug - species: clinical breakpoint. MLIS?
#   * Compare geom_smooth(method = "auto", se = TRUE) MIC tested against and evolved in with Maltas 2025 which seems to imply a downward effect but resistance to selecting drug is proxy for time and do MIC day 6 / MIC day 4, day 4 / day 2, day 2 / WT
#   * Check that MIC evolved in is actually R !! @ Maltas 2019 huh
#   * Add drug classes e.g .Barbora 2017 has annotated some, Maltas 2025 too
#         e.g. Nichol 2019: Beta–lactam antibiotics is overarching for Aminopenicillin & Cephalosporin
#   * Add clinical breakpoint (do I need MIC WT for this? Or at least MIC mutant right as concentration based)
#   * Correct Sakenova: have info on frequency of effect (not strength), but not necessary for unique drug pair percentage figure
#   * Check Durao 2015 Rifampicin and Streptomycin
#   * Repeat fig 3c of Barbosa 2017
#   * Barbosa used combination drug: PIT = PIP + ... --> duplicate to include as if only one of tro
#          = bias if comes from other drug effect, or keep as seperate entry or exclude?
# Classify by static and cidal drugs

# Questions to answer
#   * Is collateral sensitivity or cross resistance more prevalent?
#   * How well do same drug evolved / tested agree? Within species, between species?
#   * How do sizes of MIC change compare?
#   * Did they all dose as high in their experimental evolution as the final drug resistance they saw? See 2025_10_ slide 16 

# Some points:
# Effect depends on pathogen 
#   e.g. TET CIP in E. faecalis (CS or no effect, Maltas) and in E. coli (CR or no effect, Barbosa and Waller)


### Function extract repeatability, directionality ####



initial_processing_show_strains = function(dataset, MIC_list_provided = FALSE, MIC_list = c()){
  
  # MIC against drug that evolved in
  if (MIC_list_provided==FALSE){
    MIC_evolved_in = dataset %>% # evolved in same drug as tested resistance for  A-->A
      dplyr::filter(Evolved_in == Tested_against)
    
    MIC_evolved_in_list = split(MIC_evolved_in$`fold MIC change`,MIC_evolved_in$Evolved_in_strain)
  }else{
    MIC_evolved_in_list = MIC_list
  }
  
  # to find minimal MIC against drug A of any repeated evolution in drug A 
  get_min_for_pattern <- function(pat, list_obj) {
    idx <- grep(pat, names(list_obj), fixed = TRUE) # needed to disregard the strain identifier
    if (length(idx) == 0) return(NA_real_)
    min(unlist(list_obj[idx]), na.rm = TRUE)
  }
  
  # to find maximal MIC against drug A of any repeated evolution in drug A 
  get_max_for_pattern <- function(pat, list_obj) {
    idx <- grep(pat, names(list_obj), fixed = TRUE)
    if (length(idx) == 0) return(NA_real_)
    max(unlist(list_obj[idx]), na.rm = TRUE)
  }
  
  # drug A --> drug A = to what MIC did you evolve when you evolved in the drug you were tested against
  dataset = dataset %>%
    dplyr::mutate(
      MIC_evolved_in = MIC_evolved_in_list[Evolved_in_strain],
      min_MIC_when_evolved_in_tested_against = # minimal MIC of all the repeats evolved in Tested_against
        vapply(Tested_against, get_min_for_pattern, 
               numeric(1), list_obj = MIC_evolved_in_list),
      max_MIC_when_evolved_in_tested_against = 
        vapply(Tested_against, get_max_for_pattern, 
               numeric(1), list_obj = MIC_evolved_in_list)) %>%
    dplyr::group_by(Evolved_in, Tested_against) %>%
    dplyr::mutate(n_strains_evolved_in = length(unique(Evolved_in_strain))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(drug_pair = paste(Evolved_in, Tested_against, sep = '->'))
  
  # A-->B and B-->A are present?
  # Unique transitions
  unique_transitions <- unique(dataset$drug_pair)
  
  # Reverse each transition
  reverse_transitions <- sub("^(.*)->(.*)$", "\\2->\\1",
                             unique_transitions)
  
  # Find transitions whose reverse is absent
  missing_reverse <- unique_transitions[
    !reverse_transitions %in% unique_transitions
  ]
  
  print('missing reverse')
  print(missing_reverse)
  dataset_filtered = dataset %>% 
    dplyr::filter(!drug_pair %in% missing_reverse)
  
  print(paste('lost rows ', nrow(dataset), 'new ', nrow(dataset_filtered)))
  
  
  # strains evolved in if tested against
  n_strains = dataset_filtered %>%
    dplyr::group_by(Evolved_in, Tested_against) %>%
    dplyr::summarise(n_strains_evolved_in = unique(n_strains_evolved_in)) %>%#only one number, used to map n_strains_tested_against evolved in
    dplyr::mutate(pair = paste0(Evolved_in, Tested_against))
  n_strains =  split(n_strains$n_strains_evolved_in, n_strains$pair)
  
  dataset_filtered = dataset_filtered %>%
    dplyr::mutate(n_strains_tested_against = n_strains[paste0(Tested_against, Evolved_in)])
  
  dataset_filtered$n_strains_tested_against = unlist(dataset_filtered['n_strains_tested_against'])
  
  return(dataset_filtered)
}


initial_processing_do_not_show_strains = function(dataset){
  # MIC against drug that evolved in
  MIC_evolved_in = dataset %>%
    dplyr::filter(Evolved_in == Tested_against)
  
  MIC_evolved_in_list = split(MIC_evolved_in$`fold MIC change`,MIC_evolved_in$Evolved_in)
  
  get_min_for_pattern <- function(pat, list_obj) {
    idx <- grep(pat, names(list_obj), fixed = TRUE)
    if (length(idx) == 0) return(NA_real_)
    min(unlist(list_obj[idx]), na.rm = TRUE)
  }
  
  get_max_for_pattern <- function(pat, list_obj) {
    idx <- grep(pat, names(list_obj), fixed = TRUE)
    if (length(idx) == 0) return(NA_real_)
    max(unlist(list_obj[idx]), na.rm = TRUE)
  }
  
  dataset =dataset %>%
    dplyr::mutate(
      MIC_evolved_in = MIC_evolved_in_list[Evolved_in],
      min_MIC_when_evolved_in_tested_against = 
        vapply(Tested_against, get_min_for_pattern, 
               numeric(1), list_obj = MIC_evolved_in_list),
      max_MIC_when_evolved_in_tested_against = 
        vapply(Tested_against, get_max_for_pattern, 
               numeric(1), list_obj = MIC_evolved_in_list)) %>%
    dplyr::group_by(Evolved_in, Tested_against) %>%
    dplyr::mutate(n_strains_evolved_in = length(unique(Evolved_in))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(drug_pair = paste(Evolved_in, Tested_against, sep = '->'))
  
  # A-->B and B-->A are present?
  # Unique transitions
  unique_transitions <- unique(dataset$drug_pair)
  
  # Reverse each transition
  reverse_transitions <- sub("^(.*)->(.*)$", "\\2->\\1",
                             unique_transitions)
  
  # Find transitions whose reverse is absent
  missing_reverse <- unique_transitions[
    !reverse_transitions %in% unique_transitions
  ]
  
  # print('unique transitions ')
  # print(unique_transitions)
  print('missing reverse')
  print(missing_reverse)
  dataset_filtered = dataset %>% 
    dplyr::filter(!drug_pair %in% missing_reverse)
  
  print(paste('lost rows ', nrow(dataset), 'new ', nrow(dataset_filtered)))
  
  n_strains = dataset_filtered %>% 
    dplyr::group_by(Evolved_in, Tested_against) %>% 
    dplyr::summarise(n_strains_evolved_in = unique(n_strains_evolved_in)) %>%
    dplyr::mutate(pair = paste0(Evolved_in, Tested_against))
  
  n_strains =  split(n_strains$n_strains_evolved_in, n_strains$pair)
  
  dataset_filtered = dataset_filtered %>% 
    dplyr::mutate(n_strains_tested_against = n_strains[paste0(Tested_against, Evolved_in)])
  dataset_filtered$n_strains_tested_against = unlist(dataset_filtered['n_strains_tested_against'])
  
  return(dataset_filtered)
}



get_min_for_pattern <- function(pat, list_obj) {
  idx <- grep(pat, names(list_obj), fixed = TRUE) # needed to disregard the strain identifier
  if (length(idx) == 0) return(NA_real_)
  min(unlist(list_obj[idx]), na.rm = TRUE)
}

# to find maximal MIC against drug A of any repeated evolution in drug A 
get_max_for_pattern <- function(pat, list_obj) {
  idx <- grep(pat, names(list_obj), fixed = TRUE)
  if (length(idx) == 0) return(NA_real_)
  max(unlist(list_obj[idx]), na.rm = TRUE)
}

combining_R_S_MIC = function(dataset, R_threshold, S_threshold){
  
  # Cross Resistance
  data_long_R = dataset %>%
    dplyr::filter(`fold MIC change` >= R_threshold & 
                    Evolved_in != Tested_against) %>% # not tested against same drug evolved in
    dplyr::group_by(Tested_against, Evolved_in) %>%
    dplyr::mutate(n_crossresistant = n(), # per direction A->B, how many are cross-resistant, can then comapre to n_strains_evolved in to check if repeatable          
                  perc_crossR = n_crossresistant/n_strains_evolved_in,
                  repeatableR_in_this_direction = case_when(
                    perc_crossR == 1 ~ TRUE, .default=FALSE),
                  repeats_evolved_in = case_when(
                    n_strains_evolved_in>1 ~ TRUE, .default=FALSE
                  )
    ) %>%
    dplyr::ungroup()
  
  # Collateral sensitivity 
  data_long_S = dataset %>%
    dplyr::filter(`fold MIC change` <= S_threshold & 
                    Evolved_in != Tested_against) %>% # not tested against same drug evolved in
    dplyr::group_by(Tested_against, Evolved_in) %>%
    dplyr::mutate(n_collaterallysensitive = n(),             
                  perc_collS = n_collaterallysensitive/n_strains_evolved_in,
                  repeatableS_in_this_direction = case_when(
                    perc_collS == 1 ~ TRUE, .default=FALSE),
                  repeats_evolved_in = case_when(
                    n_strains_evolved_in>1 ~ TRUE, .default=FALSE)
    ) %>%
    dplyr::ungroup()
  
  data_long_notRS = dataset %>%
    dplyr::filter(`fold MIC change` < R_threshold & 
                    `fold MIC change` > S_threshold &
                    Evolved_in != Tested_against) %>% # not tested against same drug evolved in
    dplyr::group_by(Tested_against, Evolved_in) %>%
    dplyr::mutate(n_noEffect = n(),             
                  perc_noEffect = n_noEffect/n_strains_evolved_in,
                  repeatablenoEffect_in_this_direction = case_when(
                    perc_noEffect == 1 ~ TRUE, .default=FALSE),
                  repeats_evolved_in = case_when(
                    n_strains_evolved_in>1 ~ TRUE, .default=FALSE)
    ) %>%
    dplyr::ungroup()
  
  # to rbind fill with NA, later will fill up/down per drug pair directionally drug A --> drug B
  data_long_R[setdiff(colnames(data_long_S), colnames(data_long_R))]=NA
  data_long_S[setdiff(colnames(data_long_R), colnames(data_long_S))]=NA
  data_long_notRS[setdiff(colnames(data_long_S), colnames(data_long_notRS))]=NA
  data_long_R[setdiff(colnames(data_long_notRS), colnames(data_long_R))]=NA
  data_long_S[setdiff(colnames(data_long_notRS), colnames(data_long_S))]=NA
  
  data_long_RS = rbind(data_long_R, data_long_S, data_long_notRS)
  
  # sometimes cross resistant, sometimes collaterally sensitive
  data_long_RS = data_long_RS %>% 
    dplyr::group_by(Evolved_in, Tested_against) %>%
    dplyr::mutate(
      Sometimes_collaterally_S_R = case_when( 
        (any(!is.na(n_collaterallysensitive)) & any(!is.na(n_crossresistant))) ~ TRUE, 
        .default=FALSE)
    ) %>%
    # fill within group to see per drug pair over all repeats which patterns emerge
    tidyr::fill(c('n_noEffect', 'perc_noEffect', 'repeatablenoEffect_in_this_direction', 
                  'repeats_evolved_in',
                  'n_collaterallysensitive', 'perc_collS', 'repeatableS_in_this_direction', 
                  'n_crossresistant', 'perc_crossR', 'repeatableR_in_this_direction'),
                .direction = 'downup') %>%
    arrange( Evolved_in, Tested_against ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      dplyr::across(
        c(perc_noEffect, perc_collS, perc_crossR, 
          n_noEffect, n_collaterallysensitive, n_crossresistant),
        ~ tidyr::replace_na(.x, 0) # as could be that some drug pairs are missing some effects and will thus remain NA unless 0
      )
    ) # means that does not exists for this drug pair
  # is missing as define per subset R, S or no effect definitions
  
  return( data_long_RS )
}




collateral_effects_characterization = function(dataset){
  
  dataset = dataset %>%
    dplyr::mutate(reciprocalCrossR = NA, # so TRUE = reciprocal, FALSE = non-reciprocal, NA = no CR in any direction
                  reciprocalCollS = NA,
                  reciprocalNoEffect = NA,
                  n_reciprocalR = 0,
                  n_reciprocalS = 0,
                  perc_reciprocalR = 0,
                  perc_reciprocalS = 0,
                  n_reciprocalnoEffect = 0,
                  perc_reciprocalnoEffect = 0)
  iteri=1
  drugs_evolved_tested = unique(dataset$Evolved_in) # already cleaned that have both directions A-->B and B-->A in initial processing
  total_drugs = length(drugs_evolved_tested)
  
  fillAB = function(dataset, column, value){
    dataset[dataset['Evolved_in']==drugA & dataset['Tested_against']==drugB,][column] = value
    return(dataset)
  }
  
  fillBA = function(dataset, column, value){
    dataset[dataset['Evolved_in']==drugB & dataset['Tested_against']==drugA,][column] = value
    return(dataset)
  }
  
  for (drugA in drugs_evolved_tested){
    for (drugB in drugs_evolved_tested[iteri:total_drugs]){
      dirAB = dataset %>% dplyr::filter(Evolved_in == drugA & Tested_against == drugB) # evolved in A--> tested against B
      dirBA = dataset %>% dplyr::filter(Tested_against == drugA & Evolved_in == drugB) # B-->A
      
      # reciprocal CR
      # if any cross resistance in both directions
      if ( any(dirAB$n_crossresistant>0, na.rm = T) & any(dirBA$n_crossresistant>0, na.rm=T) ){
        # reciprocal status
        dataset = fillAB(dataset, 'reciprocalCrossR', TRUE)
        dataset = fillBA(dataset, 'reciprocalCrossR', TRUE)
      }
      
      # other direction CR details
      if ( any(dirAB$n_crossresistant>0, na.rm = T) ){
        # details of other direction
        # number of repeats cross R
        dataset = fillBA(dataset, 'n_reciprocalR',
                         unique(dirAB$n_crossresistant))
        # percent of repeats showing cross R in other direction
        dataset = fillBA(dataset, 'perc_reciprocalR',
                         unique(dirAB$perc_crossR))
      }
      
      if (any(dirBA$n_crossresistant>0, na.rm=T)){
        dataset = fillAB(dataset, 'n_reciprocalR',
                         unique(dirBA$n_crossresistant))
        dataset = fillAB(dataset, 'perc_reciprocalR',
                         unique(dirBA$perc_crossR))
      }
      
      # reciprocal CS
      if ( any(dirAB$n_collaterallysensitive>0, na.rm = T) & any(dirBA$n_collaterallysensitive>0, na.rm=T) ){
        dataset = fillAB(dataset, 'reciprocalCollS', TRUE)
        dataset = fillBA(dataset, 'reciprocalCollS', TRUE)
      }
      
      # other direction CS  details
      if ( any(dirAB$n_collaterallysensitive>0, na.rm = T) ){
        # details of other direction
        # number of repeats CS
        dataset = fillBA(dataset, 'n_reciprocalS',
                         unique(dirAB$n_collaterallysensitive))
        # percent of repeats showing CS in other direction
        dataset = fillBA(dataset, 'perc_reciprocalS',
                         unique(dirAB$perc_collS))
      }
      
      if (any(dirBA$n_collaterallysensitive>0, na.rm=T)){
        dataset = fillAB(dataset, 'n_reciprocalS',
                         unique(dirBA$n_collaterallysensitive))
        dataset = fillAB(dataset, 'perc_reciprocalS',
                         unique(dirBA$perc_collS))
      }
      
      
      # no effect? 
      
      if ( any(dirAB$n_noEffect>0, na.rm = T) & any(dirBA$n_noEffect>0, na.rm=T) ){
        dataset = fillAB(dataset, 'reciprocalNoEffect', TRUE)
        dataset = fillBA(dataset, 'reciprocalNoEffect', TRUE)
      }
      
      if ( any(dirAB$n_noEffect>0, na.rm=T) ){
        dataset = fillBA(dataset, 'n_reciprocalnoEffect',
                         unique(dirAB$n_noEffect))
        dataset = fillBA(dataset, 'perc_reciprocalnoEffect',
                         unique(dirAB$perc_noEffect))
      }
      if ( any(dirBA$n_noEffect>0, na.rm=T) ){
        dataset = fillAB(dataset, 'n_reciprocalnoEffect',
                         unique(dirBA$n_noEffect))
        dataset = fillAB(dataset, 'perc_reciprocalnoEffect',
                         unique(dirBA$perc_noEffect))
      }
      
      # non-reciprocal
      
      # CR in one direction not other
      if ( any(dirAB$n_crossresistant>0, na.rm = T) & !any(dirBA$n_crossresistant>0, na.rm=T)){
        # reciprocal status
        dataset = fillAB(dataset, 'reciprocalCrossR', FALSE)
        dataset = fillBA(dataset, 'reciprocalCrossR', FALSE)
      }
      
      if ( !any(dirAB$n_crossresistant>0, na.rm = T) & any(dirBA$n_crossresistant>0, na.rm=T)){
        # reciprocal status
        dataset = fillBA(dataset, 'reciprocalCrossR', FALSE)
        dataset = fillAB(dataset, 'reciprocalCrossR', FALSE)
      }
      
      # CS in one direction not other
      if ( any(dirAB$n_collaterallysensitive>0, na.rm = T) & !any(dirBA$n_collaterallysensitive>0, na.rm=T)){
        # reciprocal status
        dataset = fillAB(dataset, 'reciprocalCollS', FALSE)
        dataset = fillBA(dataset, 'reciprocalCollS', FALSE)
      }
      
      if ( !any(dirAB$n_collaterallysensitive>0, na.rm = T) & any(dirBA$n_collaterallysensitive>0, na.rm=T)){
        # reciprocal status
        dataset = fillBA(dataset, 'reciprocalCollS', FALSE)
        dataset = fillAB(dataset, 'reciprocalCollS', FALSE)
      }
      
      
      # No Effect in one direction not other
      if ( any(dirAB$n_noEffect>0, na.rm = T) & !any(dirBA$n_noEffect>0, na.rm=T)){
        # reciprocal status
        dataset = fillAB(dataset, 'reciprocalNoEffect', FALSE)
        dataset = fillBA(dataset, 'reciprocalNoEffect', FALSE)
      }
      
      if ( !any(dirAB$n_noEffect>0, na.rm = T) & any(dirBA$n_noEffect>0, na.rm=T)){
        # reciprocal status
        dataset = fillBA(dataset, 'reciprocalNoEffect', FALSE)
        dataset = fillAB(dataset, 'reciprocalNoEffect', FALSE)
      }
      
    }
    iteri=iteri+1
  }
  
  return(dataset)
  
}


### Drugs abbreviations ####



drugs_abbv = list(
  # Sakenova
  'Cefaclor' = 'CEC', 'Azithromycin' = 'AZY', 'Cefoxitin' = 'CFX',
  'Levofloxacin' = 'LEV', 'Doxycycline' = 'DOX', 'Polymyxin B' = 'POL',
  'Gentamicin'='GEN', 'Norfloxacin' = 'NOR', 'Nitrofurantoin' = 'NIT',
  'Mecillinam' = 'MEC', 'Nalidixic acid' = 'NAL', 'Sulfamethizole' = 'SUL',  # unsure about Sulfamethizole
  'Tobramycin' = 'TOB', 'Fosfomycin' = 'FOF', 'Minocycline' = 'MIN',
  'Chloramphenicol' = 'CHL', 'Ciprofloxacin' = 'CIP', 'Streptomycin' = 'STRP',
  'Trimethoprim' = 'TMP', 'Amoxicillin' = 'AMX', 'Ampicillin' = 'AMP',
  'Tetracycline' = 'TET', 'Ceftazidime' = 'CEF',  
  # Maltas 2019
  'Oxacillin' = 'OXA',
  # Maltas 2025
  'Daptomycin' = 'DAP', 'Linezolid' = 'LZD', 'Ceftriaxone' = 'CRO',
  # Barbosa
  'Gentamicin' = 'GEN', 'Streptomycin' = 'STR',
  'Piperacillin and Tazobactam' = 'PIT', 'Carbenicillin' = 'CAR',
  'Doripenem' = 'DOR', 'Imipenem' = 'IMI', 'Cefsulodin' = 'CEF'
  # Li 2026
)

# Maltas 2019 has full name already
drugs_full_sak = list('Cefaclor' = 'CEC', 'Azithromycin' = 'AZY', 'Cefoxitin' = 'CFX',
                      'Levofloxacin' = 'LEV', 'Doxycycline' = 'DOX', 'Polymyxin B' = 'POL',
                      'Gentamicin'='GEN', 'Norfloxacin' = 'NOR', 'Nitrofurantoin' = 'NIT',
                      'Mecillinam' = 'MEC', 'Nalidixic acid' = 'NAL', 'Sulfamethizole' = 'SUL',  # unsure about Sulfamethizole
                      'Tobramycin' = 'TOB', 'Fosfomycin' = 'FOF', 'Minocycline' = 'MIN',
                      'Chloramphenicol' = 'CHL', 'Ciprofloxacin' = 'CIP', 'Streptomycin' = 'STRP',
                      'Trimethoprim' = 'TMP', 'Amoxicillin' = 'AMX', 'Ampicillin' = 'AMP',
                      'Tetracycline' = 'TET', 'Ceftazidime' = 'CEF') 
drugs_full_sak = with(stack(drugs_full_sak), split(as.character(ind), values))

drugs_full_Barbosa = list('Gentamicin' = 'GEN', 'Streptomycin' = 'STR',
                          'Piperacillin and Tazobactam' = 'PIT', 'Carbenicillin' = 'CAR',
                          'Doripenem' = 'DOR', 'Imipenem' = 'IMI', 'Cefsulodin' = 'CEF',
                          'Ciprofloxacin' = 'CIP', 'Control' = 'CON') # control will be filtered out as never tested against this 
drugs_full_Barbosa = with(stack(drugs_full_Barbosa), split(as.character(ind), values))

drug_name_Waller <- c(
  "Rifampicin",
  "Fidaxomicin",
  "Linezolid",
  "Kanamycin",
  "Streptomycin",
  "Capreomycin",
  "Levofloxacin",
  "Bedaquiline",
  "Clofazimine",
  "Q203",
  "TB-47",
  "Pretomanid",
  "Isoniazid",
  "Ethionamide",
  "Ethambutol",
  "PBTZ-169",
  "AZ7371",
  "SQ109",
  "Thioacetazone",
  "Thiocarlide",
  "Thiophene-2",
  "Compound-1 (Molecule #429)",
  "Tunicamycin",
  "Nitrofurantoin"
)

abbreviation <- c(
  "RIF", "FIX", "LZD", "KAN", "STRP", "CAP",
  "LEV", "BDQ", "CFZ", "Q203", "TB47", "PA824",
  "INH", "ETH", "EMB", "PBTZ", "AZ7", "SQ109",
  "TAC", "TCL", "THP", "Cpd1", "TUN", "NFT"
)
drugs_full_Waller = setNames(as.list(drug_name_Waller), abbreviation)

drug_name_Maltas2025 <- c(
  "Ceftriaxone",
  "Ciprofloxacin",
  "Daptomycin",
  "Doxycycline",
  "Linezolid"
)

abbreviation <- c(
  "CRO",
  "CIP",
  "DAP",
  "DOX",
  "LZD"
)
drugs_full_Maltas2025 = setNames(as.list(drug_name_Maltas2025), abbreviation)

# missin Imamovic and Li

antibiotic_Imamovic2013 <- c(
  "Amikacin",
  "Gentamicin",
  "Kanamycin",
  "Streptomycin",
  "Ampicillin",
  "Amoxicillin",
  "Piperacillin",
  "Cefuroxime",
  "Cefepime",
  "Nalidixic acid",
  "Ciprofloxacin",
  "Levofloxacin",
  "Tetracycline",
  "Minocycline",
  "Tigecycline",
  "Chloramphenicol",
  "Azithromycin",
  "Colistin",
  "Polymyxin B",
  "Fosfomycin",
  "Rifampicin",
  "Nitrofurantoin",
  "Trimethoprim"
)

abbreviation <- c(
  "AMI", "GEN", "KAN", "STR",
  "AMP", "AMX", "PIP", "CFX", "CFP",
  "NAL", "CIP", "LEV",
  "TET", "MIN", "TGC",
  "CHL", "AZY", "COL", "POL",
  "FOS", "RIF", "NIT", "TRI"
)

drugs_full_Imamovic2013 = setNames(as.list(antibiotic_Imamovic2013), abbreviation)

antibiotic_Imamovic2018 <- c(
  "Amikacin",
  "Gentamicin",
  "Tobramycin",
  "Ciprofloxacin",
  "Levofloxacin",
  "Ampicillin",
  "Piperacillin",
  "Carbenicillin",
  "Ticarcillin",
  "Aztreonam",
  "Cefepime",
  "Cefuroxime",
  "Ceftazidime",
  "Meropenem",
  "Imipenem",
  "Minocycline",
  "Doxycycline",
  "Azithromycin",
  "Erythromycin",
  "Clarithromycin",
  "Colistin",
  "Fosfomycin",
  "Rifampicin",
  "Trimethoprim/Sulfamethoxazole"
)

abbreviation <- c(
  "AMI", "GEN", "TOB", "CIP", "LEV",
  "AMP", "PIP", "CAR", "TIC", "AZE",
  "CFP", "CFX", "CFZ", "MER", "IMI",
  "MIN", "DOX", "AZY", "ERI", "CLA",
  "COL", "FOS", "RIF", "TMS"
)

drugs_full_Imamovic2018 = setNames(as.list(antibiotic_Imamovic2018), abbreviation)

antibiotic_Li2026 <- c(
  "Cefotaxime", "Ceftiofur", "Meropenem", "Levofloxacin", "Ciprofloxacin",
  "Apramycin", "Amikacin", "Doxycycline", "Tetracycline", "Chloramphenicol",
  "Florfenicol", "Sulfamethoxazole/Trimethoprim", "Fosfomycin", "Colistin"
)

abbreviation <- c(
  "CTX", "CEF", "MEM", "LVF", "CIP", "APR", "AMK", "DOX", "TET", "CHL",
  "FFC", "SMX/TMP", "FOS", "COL"
)

drugs_full_Li2026 = setNames(as.list(antibiotic_Li2026), abbreviation)

drug_class_Li2026 <- c(
  "Cephalosporin", "Cephalosporin", "Carbapenem", "Quinolone", "Quinolone",
  "Aminoglycoside", "Aminoglycoside", "Tetracycline", "Tetracycline",
  "Amphenicol", "Amphenicol", "Sulfonamide", "Fosfomycin", "Polymyxin"
)

drug_name_Maltas2019 <- c(
  "Daptomycin", "Ampicillin", "Oxacillin", "Ceftriaxone", "Fosfomycin",
  "Tetracycline", "Doxycycline", "Tigecycline", "Spectinomycin",
  "Linezolid", "Chloramphenicol", "Ciprofloxacin", "Levofloxacin",
  "Nitrofurantoin", "Rifampicin"
)

drug_class_Maltas2019 <- c(
  "Lipopeptide", "Beta-Lactam", "Beta-Lactam", "Beta-Lactam", "Fosfomycin",
  "Tetracycline", "Tetracycline", "Tetracycline", "Aminoglycosides",
  "Oxazolidinone", "Amphenicol", "Quinolone", "Quinolone",
  "Nitrofuran", "Rifamycin"
)

mechanism_of_action_Maltas2019 <- c(
  "Cell membrane insertion",
  "Inhibits cell wall synthesis",
  "Inhibits cell wall synthesis",
  "Inhibits cell wall synthesis",
  "Inhibits cell wall synthesis",
  "30S protein synthesis inhibitor",
  "30S protein synthesis inhibitor",
  "30S protein synthesis inhibitor",
  "30S protein synthesis inhibitor",
  "50S protein synthesis inhibitor",
  "50S protein synthesis inhibitor",
  "DNA gyrase inhibitor",
  "DNA gyrase inhibitor",
  "Multiple mechanisms",
  "RNA polymerase inhibitor"
)

target_Li2026 <- c(
  "Cell wall",
  "Cell wall",
  "Cell wall",
  "DNA gyrase",
  "DNA gyrase",
  "Protein synthesis, 30S",
  "Protein synthesis, 30S",
  "Protein synthesis, 30S",
  "Protein synthesis, 30S",
  "Protein synthesis, 50S",
  "Protein synthesis, 50S",
  "Folic acid synthesis",
  "Cell wall biogenesis",
  "Lipopolysaccharide"
)

drug_class_Imamovic2018 <- c(
  "aminoglycoside",
  "aminoglycoside",
  "aminoglycoside",
  "quinolone",
  "quinolone",
  "b-lactam (penicillin)",
  "b-lactam (penicillin)",
  "b-lactam (penicillin)",
  "b-lactam (penicillin)",
  "b-lactam (monobactam)",
  "b-lactam (cephalosporin)",
  "b-lactam (cephalosporin)",
  "b-lactam (cephalosporin)",
  "b-lactam (carbapenem)",
  "b-lactam (carbapenem)",
  "tetracycline",
  "tetracycline",
  "macrolide",
  "macrolide",
  "macrolide",
  "polymyxin",
  "fosfomycin",
  "rifamycin",
  "antifolate"
)

target_Imamovic2018 <- c(
  "30S ribosomal subunit",
  "30S ribosomal subunit",
  "30S ribosomal subunit",
  "DNA gyrase / topoisomerase IV",
  "DNA gyrase / topoisomerase IV",
  "Cell wall synthesis",
  "Cell wall synthesis",
  "Cell wall synthesis",
  "Cell wall synthesis",
  "Cell wall synthesis",
  "Cell wall synthesis",
  "Cell wall synthesis",
  "Cell wall synthesis",
  "Cell wall synthesis",
  "Cell wall synthesis",
  "30S ribosomal subunit",
  "30S ribosomal subunit",
  "50S ribosomal subunit",
  "50S ribosomal subunit",
  "50S ribosomal subunit",
  "Cell membrane",
  "Cell wall synthesis",
  "RNA polymerase",
  "Folate synthesis"
)

drug_class_imamovic2013 <- c(
  "Aminoglycoside",
  "Aminoglycoside",
  "Aminoglycoside",
  "Aminoglycoside",
  "Beta-Lactam",
  "Beta-Lactam",
  "Beta-Lactam",
  "Beta-Lactam",
  "Beta-Lactam",
  "Quinolone",
  "Quinolone",
  "Quinolone",
  "Tetracycline",
  "Tetracycline",
  "Tetracycline",
  "Amphenicol",
  "Macrolide",
  "Polymyxin",
  "Polymyxin",
  "Fosfomycin",
  "Rifamycin",
  "Nitrofuran",
  "Dihydrofolate reductase inhibitor"
)

target_Imamovic2013 <- c(
  "Protein synthesis, 30S",
  "Protein synthesis, 30S",
  "Protein synthesis, 30S",
  "Protein synthesis, 30S",
  "Cell wall",
  "Cell wall",
  "Cell wall",
  "Cell wall",
  "Cell wall",
  "DNA gyrase",
  "DNA gyrase",
  "DNA gyrase",
  "Protein synthesis, 30S",
  "Protein synthesis, 30S",
  "Protein synthesis, 30S",
  "Protein synthesis, 50S",
  "Protein synthesis, 50S",
  "Lipopolysaccharide",
  "Lipopolysaccharide",
  "Cell wall biogenesis",
  "RNA polymerase",
  "Multiple",
  "Folic acid biosynthesis"
)

drug_class_Waller <- c(
  "Rifamycin",             # Rifampicin
  "Macrolide",             # Fidaxomicin
  "Oxazolidinone",         # Linezolid
  "Aminoglycoside",        # Kanamycin
  "Aminoglycoside",        # Streptomycin
  "Cyclic peptide",        # Capreomycin
  "Fluoroquinolone",       # Levofloxacin
  "Diarylquinoline",       # Bedaquiline
  "Riminophenazine",       # Clofazimine
  "Imidazopyridine",       # Q203
  NA,                      # TB-47
  "Nitroimidazole",        # Pretomanid
  "Isonicotinic acid derivative", # Isoniazid
  "Thioamide",             # Ethionamide
  NA,                      # Ethambutol
  "Benzothiazinone",       # PBTZ-169
  NA,                      # AZ7371
  "Diamine derivative",    # SQ109
  "Thiosemicarbazone",     # Thioacetazone
  NA,                      # Thiocarlide
  NA,                      # Thiophene-2
  NA,                      # Compound-1 (Molecule #429)
  "Nucleoside antibiotic", # Tunicamycin
  "Nitrofuran"             # Nitrofurantoin
)

Waller_mechanism <- c(
  "Transcription",
  "Transcription",
  "Translation",
  "Translation",
  "Translation",
  "Translation",
  "DNA gyrase",
  "Bioenergetics",
  "Bioenergetics",
  "Bioenergetics",
  "Bioenergetics",
  "Bioenergetics/cell wall",
  "Cell wall",
  "Cell wall",
  "Cell wall",
  "Cell wall",
  "Cell wall",
  "Cell wall",
  "Cell wall",
  "Cell wall",
  "Cell wall",
  "Cell wall",
  "Cell wall",
  "Undefined"
)

Maltas2025_drug_class <- c(
  "Beta-Lactam",
  "Quinolone",
  "Lipopeptide",
  "Tetracycline",
  "Oxazolidinone"
)

Maltas2025_mechanism <- c(
  "Cell wall synthesis inhibitor",
  "DNA gyrase inhibitor",
  "Cell membrane insertion",
  "30S protein synthesis inhibitor",
  "50S protein synthesis inhibitor"
)

functional_target_Barbosa2017 <- c(
  "DNA repair",
  "Protein synthesis",
  "Protein synthesis",
  "Cell wall synthesis",
  "Cell wall synthesis",
  "Cell wall synthesis",
  "Cell wall synthesis",
  "Cell wall synthesis"
)

drug_class_Barbosa2017 <- c(
  "Quinolones",
  "Aminoglycosides",
  "Aminoglycosides",
  "Penicillins",
  "Penicillins",
  "Carbapenems",
  "Carbapenems",
  "Cephalosporins"
)

drug_Barbosa2017 <- c(
  "Ciprofloxacin",
  "Gentamicin",
  "Streptomycin",
  "Piperacillin and Tazobactam",
  "Carbenicillin",
  "Doripenem",
  "Imipenem",
  "Cefsulodin"
)

# processing Sakenova drug classes 
Sak_fig4 = as.data.frame(
  read_excel('/Users/louagema/Documents/data/Sakenova 2025/Sakenova_2025_fig4 number of interactions.xlsx',
             sheet=1, col_names=TRUE))
all_drugs_Sak = c(Sak_fig4$drug1, Sak_fig4$drug2)
all_classes_Sak = c(Sak_fig4$`Drug category (drug1)`, Sak_fig4$`Drug category (drug2)`)

Sak_drug_class <- setNames(
  as.list(all_classes_Sak[!duplicated(all_drugs_Sak)]),
  all_drugs_Sak[!duplicated(all_drugs_Sak)]
)

Sak_drug_class[Sak_drug_class == 'nitrofurantoin'] = 'Nitrofuran'

print(unique(Sak_fig4$drug1)[! unique(Sak_fig4$drug1) %in% unique(Sak_fig4$drug2)])
print(unique(Sak_fig4$drug2)[! unique(Sak_fig4$drug2) %in% unique(Sak_fig4$drug1)])


##### All drug Classes ####### 
all_drug_classes = c(
  setNames(as.list(drug_class_Li2026), antibiotic_Li2026),
  setNames(as.list(Maltas2025_drug_class), drug_name_Maltas2025),
  setNames(as.list(drug_class_Waller), drug_name_Waller),
  # waller did not list drug class
  setNames(as.list(drug_class_imamovic2013), antibiotic_Imamovic2013),
  setNames(as.list(drug_class_Imamovic2018), antibiotic_Imamovic2018),
  setNames(as.list(drug_class_Maltas2019), drug_name_Maltas2019),
  Sak_drug_class,
  setNames(as.list(drug_class_Barbosa2017), drug_Barbosa2017)
)

all_drug_classes = all_drug_classes[order(names(all_drug_classes))]

names(all_drug_classes) <- paste0(
  toupper(substr(names(all_drug_classes), 1, 1)),
  tolower(substring(names(all_drug_classes), 2))
)

sort(unique(names(all_drug_classes)))
# "Trimethoprim/sulfamethoxazole" both orders
# Polymyxinb with and without space, same for "Nalidixic acid"   
inter_names = names(all_drug_classes)
inter_names[inter_names == "Trimethoprim/sulfamethoxazole"] = "Sulfamethoxazole/Trimethoprim"
inter_names[inter_names == "Sulfamethoxazole/trimethoprim"] = "Sulfamethoxazole/Trimethoprim"
inter_names[inter_names == "Nalidixicacid"] = "Nalidixic acid"
inter_names[inter_names == "Polymyxinb"] = "Polymyxin b"
names(all_drug_classes) = inter_names
names(all_drug_classes)

all_drug_classes <- lapply(all_drug_classes, function(x) {
  if (!is.na(x)){
    paste0(
      toupper(substr(x, 1, 1)),
      tolower(substring(x, 2))
    )
  }
  else{
    x
  }
  
})

sort(unique(unlist(all_drug_classes, use.names=FALSE)))
inter_classes = unlist(all_drug_classes, use.names=FALSE)
all_drug_classes[all_drug_classes == 'Aminoglycosides'] = 'Aminoglycoside'
all_drug_classes = lapply(all_drug_classes, function(x) {
  gsub("\\s*\\([^\\)]+\\)", "", x)
})
all_drug_classes[all_drug_classes == "B-lactam"] = 'Beta-lactam'
all_drug_classes[all_drug_classes == "Beta-lactams"] = 'Beta-lactam'
all_drug_classes[all_drug_classes == "Beta-Lactam"] = 'Beta-lactam'
all_drug_classes[all_drug_classes == "Beta-lactam"] = 'Beta-lactam'
all_drug_classes[all_drug_classes == "Carbapenems"] = 'Carbapenem'
all_drug_classes[all_drug_classes == "Fluoroquinolone"] = 'Quinolone'
all_drug_classes[all_drug_classes == "Macrolides"] = 'Macrolide'
all_drug_classes[all_drug_classes == "Polymyxin"] = 'Polymyxin b'
all_drug_classes[all_drug_classes == "Quinolones"] = 'Quinolone'
all_drug_classes[all_drug_classes == "Sulfonamides"] = 'Sulfonamide'
all_drug_classes[all_drug_classes == "Tetracyclines"] = 'Tetracycline'


# Resistance
# following Sakenova
R_threshold = 2
S_threshold = 0.5



### Waller 2023 ####



data_raw = as.data.frame(
  read_excel('/Users/louagema/Documents/data/Waller 2023/Waller 2023 cross resistance profiles.xlsx',
             sheet=1, col_names=TRUE))
data <- data_raw[,-1]
rownames(data) <- data_raw[,1]

# which drugs are tested against
print(unique(substr(rownames(data), 1, 3)))
print(length(unique(substr(rownames(data), 1, 3))))
max(data, na.rm=TRUE)
min(data, na.rm=TRUE)

# which drugs are evolved in and tested against
common_drugs_evolved_tested = intersect(
  unique(rownames(data)),
  unique(sub('\\-.*','', colnames(data))))
print(common_drugs_evolved_tested)
print(length(common_drugs_evolved_tested))

# long format
inter = data
inter$Tested_against = rownames(inter)
data_long = inter %>%
  tidyr::pivot_longer(cols = !Tested_against, 
                      names_to = c('Evolved_in', 'strain'), names_sep = '-',
                      values_to = 'average log(3) fold MIC change') %>%
  dplyr::mutate(Evolved_in_old = Evolved_in,
                Evolved_in = unlist(drugs_full_Waller[Evolved_in], use.names = FALSE),
                Tested_against_old = Tested_against,
                Tested_against = unlist(drugs_full_Waller[Tested_against], use.names = FALSE),
                'fold MIC change' = 3^`average log(3) fold MIC change`,
                'Evolved_in_strain' = paste(Evolved_in, strain, sep = '-'))%>%
  tidyr::drop_na(`fold MIC change`) 

typeof(data_long$Evolved_in)
rownames(inter)[!rownames(inter) %in% names(drugs_full_Waller)]

# filter only those that came from same isolating concentration 
# as else not repeats of same condition and cannot use for repeatability analysis
# run this if want to run only repeats later

Waller_check_MIC_WT = as.data.frame(
  read_excel('/Users/louagema/Documents/data/Waller 2023/Waller 2023 cross resistance profiles.xlsx',
             sheet='Table1', col_names=TRUE))
inter = split(Waller_check_MIC_WT$`Liquid MIC against M tuberculosis mc26206 susceptible parent strain uM`,
              Waller_check_MIC_WT$`Abbreviation`)

Waller_check_isolatingConc = as.data.frame(
  read_excel('/Users/louagema/Documents/data/Waller 2023/Waller 2023 cross resistance profiles.xlsx',
             sheet='Table2', col_names=TRUE)) %>%
  tidyr::separate(`Resistant mutant`, c('Evolved_in', 'strain'), sep = '-') %>%
  dplyr::mutate(Evolved_in_strain = paste(Evolved_in, strain, sep = '-RM'))

inter2 = split(Waller_check_isolatingConc$`Isolating concentration uM`,
               Waller_check_isolatingConc$`Evolved_in_strain`)

data_long_prefilter = data_long %>%
  dplyr::mutate(Evolved_in_strain_old = paste(Evolved_in_old, strain, sep='-'),
                MIC_WT_Evolved_in = inter[toupper(Evolved_in_old)],
                MIC_WT_Tested_against = inter[toupper(Tested_against_old)],
                isolating_conc = inter2[toupper(Evolved_in_strain_old)]
  )

data_long = data_long_prefilter %>%
  dplyr::group_by(Evolved_in, Tested_against, isolating_conc) %>%
  dplyr::mutate(n_strains_evolved_in = n()) %>%
  # dplyr::filter(n_strains_evolved_in>1) %>%
  dplyr::ungroup()

# start processing
data_long = initial_processing_show_strains(data_long)

# correcting n_strains_tested against as some concentrations tested different number of times as others
# ND in excel
# View(data_long %>% dplyr::filter(length(n_strains_tested_against)>1))

data_long$isolating_conc = unlist(data_long$isolating_conc, use.names = FALSE)

n_strains = data_long %>%
  dplyr::group_by(Evolved_in, Tested_against, isolating_conc) %>%
  dplyr::summarise(
    n_strains_evolved_in = n_distinct(Evolved_in_strain)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(evolved_in_tested_against = paste0(Evolved_in, Tested_against))

# View(n_strains %>% dplyr::filter(length(n_strains_evolved_in)>1))

n_strains2 =  split(n_strains$n_strains_evolved_in, n_strains$evolved_in_tested_against)

data_long = data_long %>%
  dplyr::group_by(Evolved_in, Tested_against, isolating_conc) %>%
  dplyr::mutate(
    n_strains_evolved_in = n_distinct(Evolved_in_strain)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(evolved_in_tested_against = paste0(Tested_against, Evolved_in),
                n_strains_tested_against = 
                  purrr::map_dbl(evolved_in_tested_against, \(x) max(n_strains2[[x]])) # maximal number of strains that adapted to drug B at one concentration
  )

# View(data_long %>% dplyr::filter(length(n_strains_tested_against)>1))

# at most 3 different concentrations of one drug that the parental strain adapted to

# maybe take pragmatic decision: keep drug concentration with most repeats
#   if equal number of repeats among drug concentrations, keep maximal drug concentration adapted to strain 

#Commented out this part for Maltas dataset
data_long_old = data_long
data_long_old %>% dplyr::group_by(n_strains_evolved_in) %>% dplyr::summarise(n=n())
data_long = data_long %>%
  dplyr::group_by(Evolved_in, Tested_against) %>%
  dplyr::filter(case_when(
    max(n_strains_evolved_in)>1 ~ n_strains_evolved_in == max(n_strains_evolved_in),
    max(n_strains_evolved_in) == 1 ~ isolating_conc == max(isolating_conc))
  )

print(paste(
  'Waller all data nrows ', nrow(data_long_old), 'now filtered to', nrow(data_long)))

print(paste(
  'Waller all data strains Evolved ', length(unique(data_long_old$Evolved_in_strain)), 'now filtered to', length(unique(data_long$Evolved_in_strain))))


data_long_RS = combining_R_S_MIC(data_long, R_threshold, S_threshold)

# A-->B and B-->A are present?
# Unique transitions
unique_transitions <- unique(data_long_RS$drug_pair)

# Reverse each transition
reverse_transitions <- sub("^(.*)->(.*)$", "\\2->\\1",
                           unique_transitions)

# Find transitions whose reverse is absent
missing_reverse <- unique_transitions[
  !reverse_transitions %in% unique_transitions
]

# print('unique transitions ')
# print(unique_transitions)
print('missing reverse')
print(missing_reverse)

# stopped here for Maltas dataset (for Waller dataset)
# 
data_long_RS = collateral_effects_characterization(
  dataset = data_long_RS)

add_metadata = function( datasetRS,
                         RReference, SSpecies, SStrain_parent,
                         SSetup, CComment, DDays_evolved,
                         GGenerations_evolved, IIncrease_every_generations,
                         IIncrease_fold_antibiotic){
  datasetRS = datasetRS %>%
    dplyr::mutate(reference = RReference,
                  Species = SSpecies,
                  Strain_parent = SStrain_parent,
                  setup = SSetup,
                  comment = CComment,
                  days_evolved =DDays_evolved,
                  generations_evolved = GGenerations_evolved,
                  increase_every_generations = IIncrease_every_generations,
                  increase_fold_antibiotic = IIncrease_fold_antibiotic)
  return(datasetRS)
}

data_long_RS = add_metadata(
  data_long_RS,
  RReference = 'Waller_2023',
  SSpecies = 'M. tuberculosis',
  SStrain_parent = 'mc26206',
  SSetup = 'Experimental evolution',
  CComment = 'Unclear what unshown repeats mean, highest concentration of compound to column 2 and diluted 3-fold along a 9-point dilution and fit but OD600',
  DDays_evolved = NA_integer_,
  GGenerations_evolved = NA_integer_,
  IIncrease_every_generations = NA_integer_,
  IIncrease_fold_antibiotic = '3point, 10-fold gradient')


# are all strains tested against all drugs??
# 7 drug pairs missed, noted with X in figure 1
not_all_waller = data_long_RS %>% 
  dplyr::group_by(Evolved_in_strain) %>% 
  dplyr::summarise(n = length(unique(Tested_against))) %>% 
  dplyr::filter(n != 19)


### Imamovic 2013 ####



data_raw_Im = as.data.frame(
  read_excel('/Users/louagema/Documents/data/Imamovic/Imamovic_2023_manual.xlsx',
             col_names=TRUE))
data_comp_Im = data_raw_Im[,2:ncol(data_raw_Im)]
rownames(data_comp_Im) = colnames(data_comp_Im)

data_automatic_Im = as.data.frame(read.csv('/Users/louagema/Documents/data/Imamovic/Data_automatic.csv'))
data_automatic_Im = data_automatic_Im[,2:ncol(data_automatic_Im)]
rownames(data_automatic_Im) = rownames(data_comp_Im)
colnames(data_automatic_Im) = colnames(data_comp_Im)

sum(data_comp_Im != data_automatic_Im)

data_Im = data_automatic_Im
data_Im = t(data_Im) %>%
  as.data.frame()

drugs_Im = colnames(data_Im)

data_Im_binary = as.matrix(data_Im)
data_Im_binary[data_Im_binary>0] = 1
data_Im_binary[data_Im_binary<0] = -1
isSymmetric(data_Im_binary) # also checks no effect so can have no effect in one direction and collateral sens or cross R in other

data_Im$Tested_against = rownames(data_Im)

data_Im_long = data_Im %>%
  tidyr::pivot_longer(cols = !Tested_against, 
                      names_to = 'Evolved_in',
                      values_to = 'fold MIC change') %>%
  dplyr::mutate(Evolved_in_old = Evolved_in,
                Evolved_in = unlist(drugs_full_Imamovic2013[Evolved_in], use.names = FALSE),
                Tested_against_old = Tested_against,
                Tested_against = unlist(drugs_full_Imamovic2013[Tested_against], use.names = FALSE),
                `fold MIC change` = case_when(
                  `fold MIC change`<0 ~ 1/abs(`fold MIC change`),
                  `fold MIC change`==0 ~ 1,
                  .default = `fold MIC change`))

typeof(data_Im_long$Evolved_in)
sum(is.na(data_Im_long$Evolved_in))
sum(is.na(data_Im_long$Tested_against))

data_Im_long = initial_processing_do_not_show_strains(data_Im_long)

data_Im_long_RS = combining_R_S_MIC(data_Im_long, R_threshold, S_threshold)

data_Im_long_RS = collateral_effects_characterization(
  dataset = data_Im_long_RS)

# A-->B and B-->A are present?
# Unique transitions
unique_transitions <- unique(data_Im_long_RS$drug_pair)

# Reverse each transition
reverse_transitions <- sub("^(.*)->(.*)$", "\\2->\\1",
                           unique_transitions)

# Find transitions whose reverse is absent
missing_reverse <- unique_transitions[
  !reverse_transitions %in% unique_transitions
]
print('missing reverse')
print(missing_reverse)


data_Im_long_RS = add_metadata(
  data_Im_long_RS,
  RReference = 'Imamovic_2013',
  SSpecies = 'E. coli',
  SStrain_parent = 'MG1655',
  SSetup = 'Experimental evolution',
  CComment = 'Do not show repeats, logarithmic drug concentration gradient with twofold dilutions, used IC90 from lowest concentration to inhibit 90%, manually extracted data from heatmap',
  DDays_evolved = 10,
  GGenerations_evolved = 100, # not explicitly mentioned but E. coli, similar medium, same temp as Sakenova (has BW25113 not MG1655 and 5 days 50 gen so maybe 100 gen))
  IIncrease_every_generations = NA_integer_,
  IIncrease_fold_antibiotic = 'Gradient plate to within transfer time can be exposed to higher concentration')


#### Imamovic clinical data EC3770 ####



data_raw_ImClin = as.data.frame(
  read_excel('/Users/louagema/Documents/data/Imamovic/Imamovic_fig6.xlsx',
             col_names=TRUE))

data_automatic_ImClin = as.data.frame(
  read.csv('/Users/louagema/Documents/data/Imamovic/Data_automatic_clin.csv')) %>%
  tidyr::drop_na()
data_automatic_ImClin = data_automatic_ImClin[,2:ncol(data_automatic_ImClin)]
rownames(data_automatic_Im) = rownames(data_comp_Im)
colnames(data_automatic_Im) = colnames(data_comp_Im)

# 4 and 8 (at least sometimes), in non clin was 4 and 2
sum(data_raw_ImClin[,2:ncol(data_raw_ImClin)] != data_automatic_ImClin)
# View(data_raw_ImClin[,2:ncol(data_raw_ImClin)] - data_automatic_ImClin)

data_ImClin = data_automatic_ImClin
rownames(data_ImClin) = data_raw_ImClin[,1]
colnames(data_ImClin) = colnames(data_raw_ImClin[,2:ncol(data_raw_ImClin)])

data_ImClin = t(data_ImClin) %>%
  as.data.frame()

drugs_ImClin = rownames(data_Im)
common_waller_imClin = intersect(drugs_ImClin, common_drugs_evolved_tested)

data_ImClin$Tested_against = rownames(data_ImClin)

data_ImClin_long = data_ImClin %>%
  tidyr::pivot_longer(cols = !Tested_against, 
                      names_to = c('strain', 'Evolved_in'), names_sep = ' ',
                      values_to = 'fold MIC change') %>%
  dplyr::filter(strain == 'EC3770') %>%
  dplyr::mutate(Evolved_in_old = Evolved_in,
                Evolved_in = unlist(drugs_full_Imamovic2013[Evolved_in], use.names = FALSE),
                Tested_against_old = Tested_against,
                Tested_against = unlist(drugs_full_Imamovic2013[Tested_against], use.names = FALSE),
                'Evolved_in_strain' = paste(Evolved_in, strain, sep = '-'),
                `fold MIC change` = case_when(
                  `fold MIC change`<0 ~ 1/abs(`fold MIC change`),
                  `fold MIC change`==0 ~ 1,
                  .default = `fold MIC change`)
  )

typeof(data_ImClin_long$Evolved_in)
sum(is.na(data_ImClin_long$Evolved_in))
sum(is.na(data_ImClin_long$Tested_against))

data_ImClin_long = initial_processing_show_strains(data_ImClin_long)

data_ImClin_long_RS = combining_R_S_MIC(data_ImClin_long, R_threshold, S_threshold)

data_ImClin_long_RS = collateral_effects_characterization(
  dataset = data_ImClin_long_RS)


# A-->B and B-->A are present?
# Unique transitions
unique_transitions <- unique(data_ImClin_long_RS$drug_pair)

# Reverse each transition
reverse_transitions <- sub("^(.*)->(.*)$", "\\2->\\1",
                           unique_transitions)

# Find transitions whose reverse is absent
missing_reverse <- unique_transitions[
  !reverse_transitions %in% unique_transitions
]
print('missing reverse')
print(missing_reverse)


data_ImClin_long_RS = add_metadata(
  data_ImClin_long_RS,
  RReference = 'Imamovic_2013',
  SSpecies = 'E. coli',
  SStrain_parent = 'EC3770',
  SSetup = 'Clinical strains evolved in vitro',                
  CComment = 'With different strains shown, logarithmic drug concentration gradient with twofold dilutions, used IC90 from lowest concentration to inhibit 90%, manually extracted data from heatmap',
  DDays_evolved = 10,
  GGenerations_evolved = 100, # not explicitly mentioned but E. coli, similar medium, same temp as Sakenova (has BW25113 not MG1655 and 5 days 50 gen so maybe 100 gen))
  IIncrease_every_generations = NA_integer_,
  IIncrease_fold_antibiotic = 'Gradient plate to within transfer time can be exposed to higher concentration')

#### Imamovic clinical data EC3856 ####



data_ImClin_longEC3856 = data_ImClin %>%
  tidyr::pivot_longer(cols = !Tested_against, 
                      names_to = c('strain', 'Evolved_in'), names_sep = ' ',
                      values_to = 'fold MIC change') %>%
  dplyr::filter(strain == 'EC3856') %>%
  dplyr::mutate(Evolved_in_old = Evolved_in,
                Evolved_in = unlist(drugs_full_Imamovic2013[Evolved_in], use.names = FALSE),
                Tested_against_old = Tested_against,
                Tested_against = unlist(drugs_full_Imamovic2013[Tested_against], use.names = FALSE),
                'Evolved_in_strain' = paste(Evolved_in, strain, sep = '-'),
                `fold MIC change` = case_when(
                  `fold MIC change`<0 ~ 1/abs(`fold MIC change`),
                  `fold MIC change`==0 ~ 1,
                  .default = `fold MIC change`)
  )

typeof(data_ImClin_longEC3856$Evolved_in)
sum(is.na(data_ImClin_longEC3856$Evolved_in))
sum(is.na(data_ImClin_longEC3856$Tested_against))


data_ImClin_longEC3856 = initial_processing_show_strains(data_ImClin_longEC3856)

data_ImClin_longEC3856_RS = combining_R_S_MIC(data_ImClin_longEC3856, R_threshold, S_threshold)

data_ImClin_longEC3856_RS = collateral_effects_characterization(
  dataset = data_ImClin_longEC3856_RS)


# A-->B and B-->A are present?
# Unique transitions
unique_transitions <- unique(data_ImClin_longEC3856_RS$drug_pair)

# Reverse each transition
reverse_transitions <- sub("^(.*)->(.*)$", "\\2->\\1",
                           unique_transitions)

# Find transitions whose reverse is absent
missing_reverse <- unique_transitions[
  !reverse_transitions %in% unique_transitions
]
print('missing reverse')
print(missing_reverse)

data_ImClin_longEC3856_RS = add_metadata(
  data_ImClin_longEC3856_RS,
  RReference = 'Imamovic_2013',
  SSpecies = 'E. coli',
  SStrain_parent = 'EC3856',
  SSetup = 'Clinical strains evolved in vitro',                
  CComment = 'With different strains shown, logarithmic drug concentration gradient with twofold dilutions, used IC90 from lowest concentration to inhibit 90%, manually extracted data from heatmap',
  DDays_evolved = 10,
  GGenerations_evolved = 100, # not explicitly mentioned but E. coli, similar medium, same temp as Sakenova (has BW25113 not MG1655 and 5 days 50 gen so maybe 100 gen))
  IIncrease_every_generations = NA_integer_,
  IIncrease_fold_antibiotic = 'Gradient plate to within transfer time can be exposed to higher concentration')


### Imamovic 2018 ####



data_raw_Im2018 = as.data.frame(
  read_excel('/Users/louagema/Documents/data/Imamovic/Imamovic_2018.xlsx',
             col_names=TRUE))

data_raw_Im2018 = data_raw_Im2018[6:32, 2:ncol(data_raw_Im2018)]
colnames(data_raw_Im2018) = data_raw_Im2018[1,]
rownames(data_raw_Im2018) = data_raw_Im2018[,1]
data_raw_Im2018 = data_raw_Im2018[2:27, 4:27] %>%
  mutate_all(as.numeric)
# fold MIC change by dividing by wt MIC
data_Im2018 = mapply('/', data_raw_Im2018[3:26,], data_raw_Im2018[1,]) %>% 
  as.data.frame()
rownames(data_Im2018) = rownames(data_raw_Im2018)[3:26]

drugs_Im2018 = colnames(data_Im2018)

data_Im2018$Evolved_in = rownames(data_Im2018)

data_Im2018_long = data_Im2018 %>%
  dplyr::mutate(Evolved_in = gsub('PAO1-','', Evolved_in))%>%
  tidyr::pivot_longer(cols = !Evolved_in, 
                      names_to = 'Tested_against',
                      values_to = 'fold MIC change') %>%
  dplyr::mutate(Evolved_in_old = Evolved_in,
                Evolved_in = unlist(drugs_full_Imamovic2018[Evolved_in], use.names = FALSE),
                Tested_against_old = Tested_against,
                Tested_against = unlist(drugs_full_Imamovic2018[Tested_against], use.names = FALSE)
  )


typeof(data_Im2018_long$Evolved_in)
sum(is.na(data_Im2018_long$Evolved_in))
sum(is.na(data_Im2018_long$Tested_against))

data_Im2018_long = initial_processing_do_not_show_strains(data_Im2018_long)

data_Im2018_long_RS = combining_R_S_MIC(data_Im2018_long, R_threshold, S_threshold)

data_Im2018_long_RS = collateral_effects_characterization(
  dataset = data_Im2018_long_RS
)


# A-->B and B-->A are present?
# Unique transitions
unique_transitions <- unique(data_Im2018_long_RS$drug_pair)

# Reverse each transition
reverse_transitions <- sub("^(.*)->(.*)$", "\\2->\\1",
                           unique_transitions)

# Find transitions whose reverse is absent
missing_reverse <- unique_transitions[
  !reverse_transitions %in% unique_transitions
]
print('missing reverse')
print(missing_reverse)

data_Im2018_long_RS = add_metadata(
  data_Im2018_long_RS,
  RReference = 'Imamovic_2018',
  SSpecies = 'P. Aeruginosa',
  SStrain_parent = 'PAO1',
  SSetup = 'Experimental evolution',
  CComment = 'No repeats, SCFM medium containing a 2-fold dilution of antibiotics',
  DDays_evolved = 10,
  GGenerations_evolved = NA_integer_, # not explicitly mentioned but E. coli, similar medium, same temp as Sakenova (has BW25113 not MG1655 and 5 days 50 gen so maybe 100 gen))
  IIncrease_every_generations = NA_integer_,
  IIncrease_fold_antibiotic = 'SCFM medium containing a 2-fold dilution of antibiotics')


### Sakenova 2025 ####


# PERSONAL CORRESPONDENCE DATA

data_raw_Sak = as.data.frame(
  read_excel('/Users/louagema/Documents/data/Sakenova 2025/Raw_Validation_results_lineages.xlsx',
             col_names=TRUE)) %>%
  dplyr::rename(Tested_against = Antibiotic,
                Evolved_in_strain = `Antibiotic-resistant strain`) 
parts <- strsplit(data_raw_Sak$Evolved_in_strain, "\\s+(?=\\d)", perl = TRUE)

data_raw_Sak$Evolved_in <- sapply(parts, `[`, 1)
data_raw_Sak$strain <- sapply(parts, `[`, 2)

common_drugs_evolved_tested_Sak = intersect(
  unique(data_raw_Sak$Evolved_in),
  unique(data_raw_Sak$Tested_against))

# never tests against drug evolved in maybe in other suppl of same paper
data_raw_Sak %>%
  dplyr::filter(Tested_against == Evolved_in)

data_Sak= data_raw_Sak %>%
  dplyr::mutate(
    `fold MIC change` =  `Fold-change`,
    MIC_evolved_in = `Fold-change in resistance`)

typeof(data_Sak$Evolved_in)
sum(is.na(data_Sak$Evolved_in))
sum(is.na(data_Sak$Tested_against))



data_Sak$Tested_against = unlist(data_Sak$Tested_against, use.names = FALSE)
data_Sak$Evolved_in = unlist(data_Sak$Evolved_in, use.names = FALSE)
data_Sak$Evolved_in_strain = paste(data_Sak$Evolved_in, data_Sak$strain, sep='-')

inti = split(data_Sak$`MIC_evolved_in`,data_Sak$Evolved_in_strain) 
inti = lapply(inti, '[[', 1)

data_Sak_long = initial_processing_show_strains(
  data_Sak, MIC_list_provided = TRUE, 
  MIC_list = inti)

data_Sak_long_RS = combining_R_S_MIC(data_Sak_long, R_threshold, S_threshold)

data_Sak_long_RS = collateral_effects_characterization(
  dataset = data_Sak_long_RS) %>%
  dplyr::mutate(MIC_evolved_in = NA_real_)


# A-->B and B-->A are present?
# Unique transitions
unique_transitions <- unique(data_Sak_long_RS$drug_pair)

# Reverse each transition
reverse_transitions <- sub("^(.*)->(.*)$", "\\2->\\1",
                           unique_transitions)

# Find transitions whose reverse is absent
missing_reverse <- unique_transitions[
  !reverse_transitions %in% unique_transitions
]
print('missing reverse')
print(missing_reverse)


data_Sak_long_RS = add_metadata(
  data_Sak_long_RS,
  RReference = 'Sakenova_2025',
  SSpecies = 'E. coli',
  SStrain_parent = 'BW25113',
  SSetup = 'Experimental evolution',
  CComment = 'Do not show all repeats (2 technical, 12 strains), eight concentrations in a twofold dilution gradient and fit but OD600, do not show raw data of evolved in drugs tested against all possible drugs evolved in only frequency',
  DDays_evolved = 5,
  GGenerations_evolved = 50,
  IIncrease_every_generations = 10,
  IIncrease_fold_antibiotic = '0.5-64x, 8 concentrations in this range then passage diluted version of population growing in highest concentration')

data_Sak_long_RS = data_Sak_long_RS %>%
  dplyr::rename(max_conc_evolved_in = MIC_evolved_in,
                min_max_conc_evolved_in_tested_against = min_MIC_when_evolved_in_tested_against,
                max_max_conc_evolved_in_tested_against = max_MIC_when_evolved_in_tested_against) # DID NOT MEASURE MIC SO IS LOWER BOUND OF CONCENTRATION THAT CAN GROW IN!


### Maltas 2019 ####


data_raw_raw_Mal = as.data.frame(
  read_excel('/Users/louagema/Documents/data/Maltas 2019.xlsx',
             sheet = 1, col_names = FALSE) )

data_raw_Mal = data_raw_raw_Mal[1:15, 1:61]
# maybe excel order c('Dapt',	'Amp',	'Ox',	'Ceft',	'Fos',	'Tet',	'Doxy',	'Tig',	'Spec',	'Line',	'Cam',	'Cip',	'Lev',	'Nitro',	'Rif')
# paper order 
rownames(data_raw_Mal) = c('Daptomycin', 'Ampicillin', 'Oxacillin', 'Ceftriaxone', 'Fosfomycin', 'Tetracycline', 'Doxycycline', 'Tigecycline', 'Spectinomycin', 'Linezolid', 'Chloramphenicol', 'Ciprofloxacin', 'Levofloxacin', 'Nitrofurantoin', 'Rifampicin')
strains = c('M1','M2','M3','M4')

data_raw_Mal = data_raw_Mal[,-1]
drugs_Mal = rownames(data_raw_Mal)
colnames(data_raw_Mal) = unlist(lapply(rownames(data_raw_Mal), 
                                       function(a) lapply(strains, function (b) paste(a, b, sep='-'))), 
                                recursive=FALSE)

# long format
inter = data_raw_Mal
inter$Tested_against = rownames(inter)
data_Mal_long = inter %>%
  tidyr::pivot_longer(cols = !Tested_against, 
                      names_to = c('Evolved_in', 'strain'), names_sep = '-',
                      values_to = 'Average log2 fold-change (drug 2)') %>%
  dplyr::mutate('Evolved_in_strain' = paste(Evolved_in, strain, sep = '-'),
                'fold MIC change' = 2^`Average log2 fold-change (drug 2)`)

data_Mal_long = initial_processing_show_strains(data_Mal_long)

data_Mal_long_RS = combining_R_S_MIC(data_Mal_long, R_threshold, S_threshold)

data_Mal_long_RS = collateral_effects_characterization(
  dataset = data_Mal_long_RS)


# A-->B and B-->A are present?
# Unique transitions
unique_transitions <- unique(data_Mal_long_RS$drug_pair)

# Reverse each transition
reverse_transitions <- sub("^(.*)->(.*)$", "\\2->\\1",
                           unique_transitions)

# Find transitions whose reverse is absent
missing_reverse <- unique_transitions[
  !reverse_transitions %in% unique_transitions
]
print('missing reverse')
print(missing_reverse)


data_Mal_long_RS = add_metadata(
  data_Mal_long_RS,
  RReference = 'Maltas_2019',
  SSpecies = 'E. faecalis',
  SStrain_parent = 'V583',
  SSetup = 'Experimental evolution',
  CComment = 'IC50, Linear dilution for dose response curve, then fit but OD600 and define collateral effect if 3-fold.',
  DDays_evolved = 8,
  GGenerations_evolved = 60,
  IIncrease_every_generations = 5,
  IIncrease_fold_antibiotic = '0.5, 2x, 4x, used highest with visible growth')



### Maltas 2025 ####



data_raw_raw_Mal25 = as.data.frame(
  read_excel('/Users/louagema/Documents/data/Maltas 2025/Maltas_2025.xlsx',
             sheet = 2, col_names = FALSE) )

data_raw_Mal25 = data_raw_raw_Mal25[3:23, 1:7]
colnames(data_raw_Mal25) = data_raw_Mal25[1,]
data_raw_Mal25 = data_raw_Mal25[-1,]
data_raw_Mal25 = data_raw_Mal25[,-2] %>%
  dplyr::rename(Evolved_in_strain = Replicate)
drugs_Mal25 = colnames(data_raw_Mal25)[-1]

# long format
inter = data_raw_Mal25
data_Mal25_long = inter %>%
  tidyr::pivot_longer(cols = !Evolved_in_strain, 
                      names_to = c('Tested_against'),
                      values_to = 'Average log2 fold-change (drug 2)') %>%
  separate(Evolved_in_strain, into = c('Evolved_in','strain'),sep = ' ') %>%
  dplyr::mutate(
    Evolved_in_old = Evolved_in,
    Evolved_in = unlist(drugs_full_Maltas2025[Evolved_in], use.names = FALSE),
    Tested_against_old = Tested_against,
    Tested_against = unlist(drugs_full_Maltas2025[Tested_against], use.names = FALSE),
    Evolved_in_strain = paste(Evolved_in, strain, sep = '-'),
    'Average log2 fold-change (drug 2)' = as.numeric(`Average log2 fold-change (drug 2)`),
    'fold MIC change' = 2^`Average log2 fold-change (drug 2)`)

typeof(data_Mal25_long$Evolved_in)
sum(is.na(data_Mal25_long$Evolved_in))
sum(is.na(data_Mal25_long$Tested_against))

data_Mal25_long = initial_processing_show_strains(data_Mal25_long)

data_Mal25_long_RS = combining_R_S_MIC(data_Mal25_long, R_threshold, S_threshold)

data_Mal25_long_RS = collateral_effects_characterization(
  dataset = data_Mal25_long_RS)



# A-->B and B-->A are present?
# Unique transitions
unique_transitions <- unique(data_Mal25_long_RS$drug_pair)

# Reverse each transition
reverse_transitions <- sub("^(.*)->(.*)$", "\\2->\\1",
                           unique_transitions)

# Find transitions whose reverse is absent
missing_reverse <- unique_transitions[
  !reverse_transitions %in% unique_transitions
]
print('missing reverse')
print(missing_reverse)

data_Mal25_long_RS = add_metadata(
  data_Mal25_long_RS,
  RReference = 'Maltas_2025',
  SSpecies = 'E. faecalis',
  SStrain_parent = 'V583',
  SSetup = 'Experimental evolution',
  CComment = 'Linear dilution for dose response curve, then fit but OD600 and define collateral effect if 3-fold.',
  DDays_evolved = 8,
  GGenerations_evolved = 60,
  IIncrease_every_generations = 5,
  IIncrease_fold_antibiotic = '0.5, 2x, 4x, used highest with visible growth')



### Barbosa 2017 ####



data_raw_Barbosa = as.data.frame(
  read.csv('/Users/louagema/Documents/data/Barbosa 2017/Barbosa_2017_Data_CollateralEffects_Figs1_2.csv')
)

unique(data_raw_Barbosa %>% dplyr::filter(Evolved_in == Tested_against) %>% dplyr::select(AUC_mean))

unique(data_raw_Barbosa %>% dplyr::filter(
  Evolved_in == Tested_against) %>% 
    dplyr::group_by(Evolved_in) %>% 
    dplyr::summarise(AUC = max(AUC_mean)) %>% dplyr::select(AUC))

data_Barbosa = data_raw_Barbosa %>%
  dplyr::mutate(
    Evolved_in_old = Evolved_in,
    Evolved_in = unlist(drugs_full_Barbosa[Evolved_in], use.names = FALSE),
    Tested_against_old = Tested_against,
    Tested_against = unlist(drugs_full_Barbosa[Tested_against], use.names = FALSE),
    Evolved_in_strain = paste(Evolved_in, Strain_id, sep = '-'))

typeof(data_Barbosa$Evolved_in)
sum(is.na(data_Barbosa$Evolved_in))
sum(is.na(data_Barbosa$Tested_against))

data_Barbosa_long = data_Barbosa %>% 
  dplyr::mutate('AUC difference' = AUC_mean)

common_drugs_evolved_tested_Barbosa = intersect(
  unique(data_Barbosa_long$Evolved_in),
  unique(data_Barbosa_long$Tested_against))
print(common_drugs_evolved_tested_Barbosa)

# AUC diff against drug that evolved in
AUC_evolved_in_Barbosa = data_Barbosa_long %>%
  dplyr::filter(Evolved_in == Tested_against)

AUC_evolved_in_list_Barbosa = split(AUC_evolved_in_Barbosa$`AUC difference`,AUC_evolved_in_Barbosa$Evolved_in_strain)

data_Barbosa_long =data_Barbosa_long %>%
  # dplyr::filter(Tested_against %in% common_drugs_evolved_tested_Barbosa & Evolved_in %in% common_drugs_evolved_tested_Barbosa)%>%
  dplyr::mutate(
    AUCdiff_evolved_in = AUC_evolved_in_list_Barbosa[Evolved_in_strain],
    min_AUCdiff_when_evolved_in_tested_against = 
      vapply(Tested_against, get_min_for_pattern, 
             numeric(1), list_obj = AUC_evolved_in_list_Barbosa),
    max_AUCdiff_when_evolved_in_tested_against = 
      vapply(Tested_against, get_max_for_pattern, 
             numeric(1), list_obj = AUC_evolved_in_list_Barbosa)) %>%
  dplyr::group_by(Evolved_in) %>%
  dplyr::mutate(n_strains_evolved_in = length(unique(Evolved_in_strain))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(drug_pair = paste(Evolved_in, Tested_against, sep = '->'))

# A-->B and B-->A are present?
# Unique transitions
unique_transitions <- unique(data_Barbosa_long$drug_pair)

# Reverse each transition
reverse_transitions <- sub("^(.*)->(.*)$", "\\2->\\1",
                           unique_transitions)

# Find transitions whose reverse is absent
missing_reverse <- unique_transitions[
  !reverse_transitions %in% unique_transitions
]

# print('unique transitions ')
# print(unique_transitions)
print('missing reverse')
print(missing_reverse)

data_Barbosa_long_filtered = data_Barbosa_long %>% 
  dplyr::filter(!drug_pair %in% missing_reverse)

print(paste('lost rows ', nrow(data_Barbosa_long), 'new ', nrow(data_Barbosa_long_filtered)))


# strains evolved in if tested against
n_strains = data_Barbosa_long_filtered %>%
  dplyr::group_by(Evolved_in, Tested_against) %>%
  dplyr::summarise(n_strains_evolved_in = unique(n_strains_evolved_in)) %>%#only one number, used to map n_strains_tested_against evolved in
  dplyr::mutate(pair = paste0(Evolved_in, Tested_against))
n_strains =  split(n_strains$n_strains_evolved_in,n_strains$pair)

data_Barbosa_long_filtered = data_Barbosa_long_filtered %>%
  dplyr::mutate(n_strains_tested_against = n_strains[paste0(Tested_against, Evolved_in)])

data_Barbosa_long_filtered$n_strains_tested_against = unlist(data_Barbosa_long_filtered['n_strains_tested_against'])

data_Barbosa_long = data_Barbosa_long_filtered
# Resistance
# discrete -32 -16 -8 -4 -2 90 2 4 8 16 32
R_threshold_AUC_diff = 0.5
S_threshold_AUC_diff = -0.5 # they seem to consider it this way
data_Barbosa_long_R = data_Barbosa_long %>%
  dplyr::filter(`AUC difference` > R_threshold_AUC_diff & 
                  Evolved_in != Tested_against) %>% # not tested against same drug evolved in
  dplyr::group_by(Tested_against, Evolved_in) %>%
  # unnecessary as no repeats anyway
  dplyr::mutate(n_crossresistant = n(),             
                perc_crossR = n_crossresistant/n_strains_evolved_in,
                repeatableR_in_this_direction = case_when(
                  perc_crossR == 1 ~ TRUE, .default=FALSE),
                repeats_evolved_in = case_when(
                  n_strains_evolved_in>1 ~ TRUE, .default=FALSE
                )
  )

data_Barbosa_long_S = data_Barbosa_long %>%
  dplyr::filter(`AUC difference` < S_threshold_AUC_diff & 
                  Evolved_in != Tested_against) %>% # not tested against same drug evolved in
  dplyr::group_by(Tested_against, Evolved_in) %>%
  dplyr::mutate(n_collaterallysensitive = n(),             
                perc_collS = n_collaterallysensitive/n_strains_evolved_in,
                repeatableS_in_this_direction = case_when(
                  perc_collS == 1 ~ TRUE, .default=FALSE),
                repeats_evolved_in = case_when(
                  n_strains_evolved_in>1 ~ TRUE, .default=FALSE)
  )

data_Barbosa_long_notRS = data_Barbosa_long %>%
  dplyr::filter(`AUC difference` < R_threshold_AUC_diff & 
                  `AUC difference` > S_threshold_AUC_diff &
                  Evolved_in != Tested_against) %>% # not tested against same drug evolved in
  dplyr::group_by(Tested_against, Evolved_in) %>%
  dplyr::mutate(n_noEffect = n(),             
                perc_noEffect = n_noEffect/n_strains_evolved_in,
                repeatablenoEffect_in_this_direction = case_when(
                  perc_noEffect == 1 ~ TRUE, .default=FALSE),
                repeats_evolved_in = case_when(
                  n_strains_evolved_in>1 ~ TRUE, .default=FALSE)
  ) %>%
  dplyr::ungroup()

data_Barbosa_long_R[setdiff(colnames(data_Barbosa_long_S), colnames(data_Barbosa_long_R))]=NA
data_Barbosa_long_S[setdiff(colnames(data_Barbosa_long_R), colnames(data_Barbosa_long_S))]=NA
data_Barbosa_long_notRS[setdiff(colnames(data_Barbosa_long_S), colnames(data_Barbosa_long_notRS))]=NA
data_Barbosa_long_R[setdiff(colnames(data_Barbosa_long_notRS), colnames(data_Barbosa_long_R))]=NA
data_Barbosa_long_S[setdiff(colnames(data_Barbosa_long_notRS), colnames(data_Barbosa_long_S))]=NA

data_Barbosa_long_RS = rbind(data_Barbosa_long_R, data_Barbosa_long_S, data_Barbosa_long_notRS)

# sometimes cross resistant, sometimes collaterally sensitive
data_Barbosa_long_RS = data_Barbosa_long_RS %>%
  dplyr::group_by(Evolved_in, Tested_against) %>%
  dplyr::mutate(
    Sometimes_collaterally_S_R = case_when(
      (any(!is.na(n_collaterallysensitive)) & any(!is.na(n_crossresistant))) ~ TRUE,
      .default=FALSE)
  )  %>%
  tidyr::fill(c('n_noEffect', 'perc_noEffect', 'repeatablenoEffect_in_this_direction', 'repeats_evolved_in',
                'n_collaterallysensitive', 'perc_collS', 'repeatableS_in_this_direction', 
                'n_crossresistant', 'perc_crossR', 'repeatableR_in_this_direction'),
              .direction = 'downup') %>%
  arrange( Evolved_in, Tested_against ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    dplyr::across(
      c(perc_noEffect, perc_collS, perc_crossR),
      ~ tidyr::replace_na(.x, 0)
    )
  ) # means that does not exists for this drug pair
# is missing as define per subset R, S or no effect definitions

drugs_Barbosa = unique(data_Barbosa_long_RS$Tested_against)

data_Barbosa_long_RS = collateral_effects_characterization(
  dataset = data_Barbosa_long_RS)

data_Barbosa_long_RS = add_metadata(
  data_Barbosa_long_RS,
  RReference = 'Barbosa_2017',
  SSpecies = 'P. Aeruginosa',
  SStrain_parent = 'PA14',
  SSetup = 'Experimental evolution',
  CComment = 'Growth curve AUC difference.',
  DDays_evolved = 24,
  GGenerations_evolved = 48,
  IIncrease_every_generations = 4,
  IIncrease_fold_antibiotic = '1.5')


### Li 2026 ####

# To do, homogenize drug. names

data_raw_Li2026 = as.data.frame(
  read_excel('/Users/louagema/Documents/data/Li 2026/Li_2026.xlsx',sheet = 'evolved MIC')
)

data_WTMIC_Li2026 = as.data.frame(
  read_excel('/Users/louagema/Documents/data/Li 2026/Li_2026.xlsx', sheet = 'parental_MIC')
)

colnames(data_WTMIC_Li2026) = c('Tested_against', 'WT', 'C600')
data_WTMIC_Li2026_long = data_WTMIC_Li2026 %>%
  tidyr::pivot_longer(-Tested_against, names_to = 'parent', values_to = 'raw_parental_MIC')


strain_list = c('WT', 'C600')
drug_list = c('MEM', 'CIP', 'FOS', 'DOX', "APR", 'COL')
rep_list = c('1', '2', '3')

x <- expand.grid(C = rep_list , B = strain_list, A = drug_list)
x = paste(x$B, x$C, x$A)
colnames(data_raw_Li2026)= c('Tested_against', x)

data_Li2026_long = data_raw_Li2026 %>%
  tidyr::pivot_longer(-Tested_against, names_to = 'parent_strain_drug', values_to = 'raw_evolved_MIC') %>%
  tidyr::separate(parent_strain_drug, into = c('parent', 'strain', 'Evolved_in'), sep = ' ') %>%
  merge(data_WTMIC_Li2026_long, by = c('parent', 'Tested_against')) %>%
  dplyr::mutate(
    Evolved_in_old = Evolved_in,
    Evolved_in = unlist(drugs_full_Li2026[Evolved_in], use.names = FALSE),
    Tested_against_old = Tested_against,
    Tested_against = unlist(drugs_full_Li2026[Tested_against], use.names = FALSE),
    `fold MIC change` = raw_evolved_MIC/raw_parental_MIC,
    Evolved_in_strain = paste(Evolved_in, strain, sep = '-'))

typeof(data_Li2026_long$Evolved_in)
sum(is.na(data_Li2026_long$Evolved_in))
sum(is.na(data_Li2026_long$Tested_against))


parental_strains = c('WT', 'C600')
i=0
for (sp in parental_strains){
  # start processing
  inter = data_Li2026_long %>% dplyr::filter(parent == sp)
  
  common_drugs_evolved_tested_inter = intersect(
    unique(inter$Evolved_in),
    unique(inter$Tested_against))
  
  inter = inter %>%
    dplyr::filter(Evolved_in %in% common_drugs_evolved_tested_inter & 
                    Tested_against %in% common_drugs_evolved_tested_inter)
  
  inter_long = initial_processing_show_strains(inter)
  
  inter_long_RS = combining_R_S_MIC(inter_long, R_threshold, S_threshold)
  
  inter_long_RS = collateral_effects_characterization(
    dataset = inter_long_RS)
  
  inter_long_RS = add_metadata(
    inter_long_RS,
    RReference = 'Li_2026',
    SSpecies = 'E. coli',
    SStrain_parent = sp,
    SSetup = 'Experimental evolution',
    CComment = 'Every day transfer new medium, five concentrations (1/4 to 4 × MIC) evolved mutants were subcultured for three passages in antibiotic-free LB to eliminate transient tolerance and confirm stable resistance',
    DDays_evolved = 14, 
    GGenerations_evolved = NA_integer_, 
    IIncrease_every_generations = NA_integer_,
    IIncrease_fold_antibiotic = NA_character_)
  
  if (i==0){
    data_Li26_long_RS = inter_long_RS
  }
  else{
    data_Li26_long_RS = rbind(data_Li26_long_RS, inter_long_RS)
  }
  i=i+1
}



### combine data ####



all_data = bind_rows(
  data_Im_long_RS, 
  data_ImClin_long_RS, 
  data_long_RS, 
  data_Mal_long_RS, 
  data_Mal25_long_RS, 
  data_ImClin_longEC3856_RS,
  data_Im2018_long_RS,
  data_Li26_long_RS
) 

all_data$MIC_evolved_in = unlist(all_data$MIC_evolved_in, use.names = FALSE)
all_data = bind_rows(all_data, data_Sak_long_RS)


all_data$MIC_evolved_in = as.numeric(unlist(all_data['MIC_evolved_in']))

# don't have MIC evolved in so do seperate
# TO DO UNLIST NOT WORKING FOR CLASS SOMEHOW class_evolved_in CHECK View(all_data[length(all_data$class_evolved_in)>1])
all_data = bind_rows(all_data, data_Barbosa_long_RS, data_Sak_long_RS) %>%
  dplyr::mutate(Evolved_in_strain = case_when(
    is.na(Evolved_in_strain) ~ Evolved_in,
    .default = Evolved_in_strain))%>%
  dplyr::mutate(
    drug_pair = paste(Evolved_in, Tested_against, sep = '->'),
    drug_pair_nodir = pmap_chr(list(Evolved_in,Tested_against), ~paste(sort(c(...)), collapse = "-"))
  )


### MIC fold change A-->B vs B-->B #####

inter2 = all_data %>%
  dplyr::filter(Tested_against != Evolved_in) %>% 
  tidyr::replace_na(list('n_collaterallysensitive'=0, 
                         'n_crossresistant'=0, 
                         'n_noEffect'=0,
                         'n_reciprocalR'=0,
                         'n_reciprocalS'=0,
                         'perc_reciprocalnoEffect'=0,
                         'n_reciprocalnoEffect'=0,
                         'perc_crossR'=0,
                         'perc_collS'=0)) %>%
  dplyr::mutate(
    Repeatable = factor(case_when(
      n_crossresistant>1 & perc_crossR==1 ~ 'Yes', 
      n_collaterallysensitive>1 & perc_collS==1 ~ 'Yes',
      n_noEffect>1 & perc_noEffect==1 ~ 'Yes',
      n_strains_evolved_in>1 & perc_crossR!=1 & perc_crossR>0 ~ 'No', 
      n_strains_evolved_in>1 & perc_collS!=1 & perc_collS>0 ~ 'No',
      n_strains_evolved_in == 1 ~ 'No repeats'
    ), 
    ,
    levels = c('Yes', 'No', 'No repeats')),
    Reciprocal = factor(case_when(
      n_crossresistant>0 & n_reciprocalR>0 & ( `fold MIC change`>=R_threshold |  `AUC difference`>=R_threshold_AUC_diff) ~ 'Yes',
      n_collaterallysensitive>0 & n_reciprocalS>0 & ( `fold MIC change`<=S_threshold |  `AUC difference`<=S_threshold_AUC_diff)~ 'Yes',
      n_crossresistant>0 & n_reciprocalR==0 ~ 'No',
      n_collaterallysensitive>0 & n_reciprocalS==0 ~ 'No',
      perc_noEffect==1 & perc_reciprocalnoEffect==1 ~ 'No collateral effect',
      (`fold MIC change`<R_threshold & `fold MIC change`>S_threshold) |  
        (`AUC difference`<R_threshold_AUC_diff & `AUC difference`>S_threshold_AUC_diff) ~
        'No collateral effect'
    ), levels = c('Yes', 'No', 'Mixed effect', 'No collateral effect', '?'))
  )


ggplot(inter2 %>% 
         dplyr::filter(!is.na(max_MIC_when_evolved_in_tested_against)) %>%
         dplyr::filter(`fold MIC change`>=R_threshold | 
                         `fold MIC change`<=S_threshold), 
       aes(x=max_MIC_when_evolved_in_tested_against, 
           y = `fold MIC change`,
           color= Reciprocal, shape = Repeatable)) + 
  # facet_grid(.~Reciprocal)+
  theme_bw(base_size = 16)+
  labs(#itle = 'MIC or IC90 fold change relative to parental strain',
    x= TeX('B $\\rightarrow$ B (max)'), #'Evolved in B, resistance to B (max of all repeats)',
    y = TeX('A $\\rightarrow$ B'), #'Evolved in A, collateral effect to B',
    color = 'Bidirectional', shape = 'Repeatable \n in this direction')+
  geom_hline(yintercept = 1)  + 
  # geom_hline(yintercept = R_threshold, alpha=0.5)  + 
  # geom_hline(yintercept = S_threshold, alpha=0.5)  + 
  # geom_vline(xintercept = R_threshold, alpha=0.5)+ 
  geom_point() +
  # geom_smooth(method = "auto", se = TRUE)+
  scale_x_log10(breaks = c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000)) + 
  scale_y_log10(labels = function(x) ifelse(x == 0, "0", x),
                breaks = c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000))+
  geom_abline(intercept=0, color ='darkgrey')+
  geom_abline(slope = -1, intercept=0, color ='darkgrey')

ggsave('/Users/louagema/Documents/Figures/CollaterEffectsConcept/AtoBvsBtoB.png',
       dpi = 300, height = 5, width = 8.5)


##### AUC Barbosa ######

inter2 = data_Barbosa_long_RS %>%
  dplyr::filter(Tested_against != Evolved_in) %>% 
  tidyr::replace_na(list('n_collaterallysensitive'=0, 
                         'n_crossresistant'=0, 
                         'n_noEffect'=0,
                         'n_reciprocalR'=0,
                         'n_reciprocalS'=0,
                         'perc_reciprocalnoEffect'=0,
                         'n_reciprocalnoEffect'=0,
                         'perc_crossR'=0,
                         'perc_collS'=0)) %>%
  dplyr::mutate(
    Repeatable = factor(case_when(
      n_crossresistant>1 & perc_crossR==1 ~ 'Yes', 
      n_collaterallysensitive>1 & perc_collS==1 ~ 'Yes',
      n_noEffect>1 & perc_noEffect==1 ~ 'Yes',
      n_strains_evolved_in>1 & perc_crossR!=1 & perc_crossR>0 ~ 'No', 
      n_strains_evolved_in>1 & perc_collS!=1 & perc_collS>0 ~ 'No',
      n_strains_evolved_in == 1 ~ 'No repeats'
    ), 
    ,
    levels = c('Yes', 'No', 'No repeats')),
    Reciprocal = factor(case_when(
      n_crossresistant>0 & n_reciprocalR>0 & (`AUC difference`>=R_threshold_AUC_diff) ~ 'Yes',
      n_collaterallysensitive>0 & n_reciprocalS>0 & (`AUC difference`<=S_threshold_AUC_diff)~ 'Yes',
      n_crossresistant>0 & n_reciprocalR==0 ~ 'No',
      n_collaterallysensitive>0 & n_reciprocalS==0 ~ 'No',
      perc_noEffect==1 & perc_reciprocalnoEffect==1 ~ 'No collateral effect',
      (`AUC difference`<R_threshold_AUC_diff & `AUC difference`>S_threshold_AUC_diff) ~
        'No collateral effect'
    ), levels = c('Yes', 'No', 'Mixed effect', 'No collateral effect', '?'))
  )


ggplot(inter2 %>% 
         dplyr::filter(`AUC difference`>=R_threshold_AUC_diff | 
                         `AUC difference`<=S_threshold_AUC_diff), 
       aes(x=max_AUCdiff_when_evolved_in_tested_against, 
           y = `AUC difference`,
           color= Reciprocal, shape = Repeatable)) + 
  # facet_grid(.~Reciprocal)+
  theme_bw(base_size = 16)+
  labs(#itle = 'MIC or IC90 fold change relative to parental strain',
    x= TeX('B $\\rightarrow$ B (max)'), #'Evolved in B, resistance to B (max of all repeats)',
    y = TeX('A $\\rightarrow$ B'), #'Evolved in A, collateral effect to B',
    color = 'Bidirectional', shape = 'Repeatable \n in this direction')+
  geom_hline(yintercept = 0)  + 
  geom_point() +
  geom_abline(intercept=0, color ='darkgrey')+
  geom_abline(slope = -1, intercept=0, color ='darkgrey')

ggsave('/Users/louagema/Documents/Figures/CollaterEffectsConcept/AtoBvsBtoB_Barbosa.png',
       dpi = 300, height = 5, width = 8.5)



### Mixed collateral effects ####

inter = all_data %>% 
  dplyr::filter(n_strains_evolved_in>1) %>% # e.g. Imamovic does not have repeats so can never have A-->B show sometimes CS and sometmeis CR
  tidyr::replace_na(list('n_collaterallysensitive'=0, 
                         'n_crossresistant'=0, 
                         'n_noEffect'=0))

# in one direction, how often see CS and CR vs only CS
a = inter %>%
  dplyr::group_by(reference, Species, Strain_parent) %>%
  dplyr::summarise(#n_drug_pairs_nodir=n_distinct(drug_pair_nodir),
    n_drug_pairs=n_distinct(drug_pair)) %>%
  
  # CS and sometimes CR
  
  merge(
    inter %>%
      dplyr::filter(n_collaterallysensitive > 0 & n_crossresistant>0) %>%
      dplyr::group_by(reference, Species, Strain_parent) %>%
      dplyr::summarise(n_CS_or_CR=n_distinct(drug_pair)),
    by = c('reference', 'Species', 'Strain_parent'), all.y=TRUE, all.x=TRUE
  ) %>%
  
  # CS and sometimes CR, never no effect
  
  merge(
    inter %>%
      dplyr::filter(n_collaterallysensitive > 0 & n_crossresistant>0 & n_noEffect==0) %>%
      dplyr::group_by(reference, Species, Strain_parent) %>%
      dplyr::summarise(n_CS_or_CR_nevernoCE=n_distinct(drug_pair)),
    by = c('reference', 'Species', 'Strain_parent'), all.y=TRUE, all.x=TRUE
  ) %>%
  
  # non-repeatable CS, sometimes see no effect, never CR
  
  merge(
    inter %>% 
      dplyr::filter(n_collaterallysensitive > 0 & n_noEffect>0 & n_crossresistant==0) %>%
      dplyr::group_by(reference, Species, Strain_parent) %>%
      dplyr::summarise(n_CS_noEffect=n_distinct(drug_pair)),
    by = c('reference', 'Species', 'Strain_parent'), all.y=TRUE, all.x=TRUE
  ) %>%
  
  # non-repeatable CR, sometimes see no effect, never CS
  
  merge(
    inter %>% 
      dplyr::filter(n_crossresistant > 0 & n_noEffect>0 & n_collaterallysensitive==0) %>%
      dplyr::group_by(reference, Species, Strain_parent) %>%
      dplyr::summarise(n_CR_noEffect=n_distinct(drug_pair)),
    by = c('reference', 'Species', 'Strain_parent'), all.y=TRUE, all.x=TRUE
  ) %>%
  
  # repeatable CS
  
  merge(
    inter %>% 
      dplyr::filter(n_collaterallysensitive > 0 & perc_collS==1) %>%
      dplyr::group_by(reference, Species, Strain_parent) %>%
      dplyr::summarise(n_CS_repeatable=n_distinct(drug_pair)),
    by = c('reference', 'Species', 'Strain_parent'), all.y=TRUE, all.x=TRUE
  ) %>%
  
  # repeatable CR
  
  merge(
    inter %>% 
      dplyr::filter(n_crossresistant > 0 & perc_crossR==1) %>%
      dplyr::group_by(reference, Species, Strain_parent) %>%
      dplyr::summarise(n_CR_repeatable=n_distinct(drug_pair)),
    by = c('reference', 'Species', 'Strain_parent'), all.y=TRUE, all.x=TRUE
  ) 

print(a)

