library(sf)
library(tidyverse)
library(stringi)
library("readxl")
library(ggpubr)
library(scales)

#######################################################################################
###################################### Figure 1a ######################################
#######################################################################################


# Replace 'path_to_shapefile' with the actual path to your shapefile

shape_data <- st_read("path_to_shapefile")

west_africa_countries <- c("BEN", "BKF", "CVI", "IVC", 
                           "GAM", "GHA", "GUI", "GBS", 
                           "LIR", "MLI", "MAU", "NER", 
                           "NIR", "SEN", "SIL", "TOG")

west_africa_shape_data <- shape_data %>%
  filter(UNDP %in% west_africa_countries)

#Replace 'path_to_sequencer_location' with the actual path to your sequencer location file
sequencer_locations_sf <- st_as_sf(read_excel("path_to_sequencer_location"), coords = c("longitude", "latitude"), crs = 4326)


#Replace 'Country_provided_Seqs' with the actual path to your countries that provided local sequence data
seqs_countries_sf <- st_as_sf(read_excel("Country_provided_Seqs"), coords = c("longitude", "latitude"), crs = 4326)


West_africa_map <- ggplot(data = west_africa_shape_data) +
  geom_sf(fill = "gray80", color = "black") +
  theme_void()+ 
  geom_sf(data = sequencer_locations_sf, size = 5, shape =21, fill = "blue") + 
  geom_sf(data= seqs_countries_sf, size = 3, shape = 13, fill = "lightgreen")+
  theme(legend.position = "none")

West_africa_map

#####################################################################################

#Replace 'metadata' with the actual path of your merged metadata file

waf_data <- read.csv("metadata", sep = "\t", header = T)


#######################################################################################
########### Figure 1b : Sex distribution in the dataset by country ####################
#######################################################################################

freq_gender_by_country <- waf_data %>%
  select(Country, Gender) %>%
  group_by(Country, Gender) %>%
  summarise(Count = n()) %>%
  mutate(freq_sex_percent = Count/sum(Count)*100)

ggplot(data = freq_gender_by_country) +
  geom_bar(aes(x= Gender, y= freq_sex_percent, fill = Country), stat = "identity", position = "stack") +
  facet_wrap(vars(Country)) +
  theme_classic() + 
  labs(x="", y = "Frequency of Sex (%)")+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = .5, hjust = 0.3, size = 8))


#######################################################################################
############### Figure 1c : sequences found in GISAID by country ######################
#######################################################################################


seq_per_country <- waf_data %>%
  group_by(Country) %>%
  summarise(Count = n(), Percent = n()/28763)

ggplot(data = seq_per_country) + 
  geom_bar(aes(x=Country, y=Percent, fill = Country), stat = "identity", position = "stack") + labs(x= NULL, y = "% of sequence per country") + 
  theme_classic() + ylim(0, 0.3) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 8))


#######################################################################################
################ Figure 1d : WHO variant repartition per country ######################
#######################################################################################

clade_per_country <- waf_data %>%
  group_by(Country, Clade, Variants_ed, Variant_Nextclade) %>%
  summarise(Count = n())

clade_per_country_filtered <- clade_per_country %>%
  filter(Variants_ed != "Unassinged") %>%
  filter(Variants_ed != "Unclassified")

ggplot(data = clade_per_country_filtered) +
  geom_bar(aes(x=Variants_ed, y=Count, fill = Country), stat = "identity", position = "stack") + ylim(0, 2500) + labs(x="", y= "Frequency of variants") +
  theme_classic() +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8))

#######################################################################################
##################### Figure 1e : Rate of variants in West Africa #####################
#######################################################################################

Variant <- waf_data %>%
  select(Country, Variants_ed) %>%
  filter(Variants_ed != "Unassinged") %>%
  filter(Variants_ed != "Unclassified") %>%
  group_by(Variants_ed) %>%
  summarise(Count = n(), freq = round((n()/28763*100),2))

ggplot(data = Variant) +
  geom_bar(aes(x=Variants_ed, y= freq, fill = Variants_ed), stat = "identity", position = "stack") + 
  geom_text(aes(x=Variants_ed, y= freq, label = freq), hjust = 0.5, vjust = -0.5) +
  labs(x= "", y="Frequence (%)") +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 18), axis.text.y = element_text(size = 8)) 

#######################################################################################
############## Figure 1f : Repartition of nextclade variant per country################
#######################################################################################

ggplot(data = clade_per_country_filtered) +
  geom_bar(aes(x=Variant_Nextclade, y=Count, fill = Country), stat = "identity", position = "stack") + ylim(0, 4000) + labs(x="", y= "Frequency of Clade") +
  theme_classic() + 
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 8))


#######################################################################################
######################### Figure 2a : Dynamics of Variants ############################
#######################################################################################

waf_data$Collection.date <- as.Date(waf_data$Collection.date)
waf_data$Collection.date.month <- format(waf_data$Collection.date, "%Y-%m")


clade_variant_time_filtered <- waf_data %>%
  group_by(Collection.date.month, Country, Variants_ed) %>%
  filter(Variants_ed != "Unassinged") %>%
  filter(Variants_ed != "Unclassified") %>%
  summarise(Count = n())

# Remove rows with missing or non-finite values in Collection.date.month
clade_variant_time_filtered <- clade_variant_time_filtered[complete.cases(clade_variant_time_filtered$Collection.date.month), ]


# Remove rows with missing values in Collection.date.month
clade_variant_time_filtered <- na.omit(clade_variant_time_filtered)

ggplot(clade_variant_time_filtered, aes(x = Collection.date.month, fill = Variants_ed)) +
  geom_bar(position = "fill") +
  labs(x = "", y = "Percentage") +
  theme_classic() +
  scale_x_discrete(labels = date_format("%Y-%m"), breaks = unique(clade_variant_time_filtered$Collection.date.month)) +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6))

#######################################################################################
###################### Figure 2b : Dynamics variant per country #######################
#######################################################################################

ggplot(clade_variant_time_filtered, aes(x= Collection.date.month, fill = Variants_ed)) + 
  geom_bar(position = "fill") +
  labs(x = "", y = "Percentage") +
  theme_classic() + 
  scale_x_discrete(labels = date_format("%Y-%m"), breaks = unique(clade_variant_time_filtered$Collection.date.month)) +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4), axis.text.y = element_text(size = 6)) +
  facet_wrap(vars(Country),nrow = 4)

  
#######################################################################################
### Figure 2c : Dynamics of variant per country according nextclade classification ####
#######################################################################################

Nextclade_time <- waf_data %>%
#  select(Collection.date.month, Country, Variant_Nextclade) %>%
  group_by(Collection.date.month, Country, Variant_Nextclade) %>%
  filter(Variant_Nextclade != "Unassigned") %>%
  filter(Variant_Nextclade != "Unclassified") %>%
  summarise(Count = n())

Nextclade_time <- Nextclade_time[complete.cases(Nextclade_time$Collection.date.month), ]



ggplot(Nextclade_time, aes(x= Collection.date.month, fill = Variant_Nextclade)) + 
  geom_bar(position = "fill") +
  labs(x = "", y = "Percentage") +
  theme_classic() + 
  scale_x_discrete(labels = date_format("%Y-%m"), breaks = unique(Nextclade_time$Collection.date.month)) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = .5, size = 5), 
        axis.text.y = element_text(size = 5)) +
  facet_wrap(vars(Country), nrow =4)


#######################################################################################
################# Figure 2d : Mutation distribution per SARS-CoV-2 genes ##############
#######################################################################################


waf_data_sub_clean <- waf_data %>%
  select(Virus.name, Country, AA.Substitutions)

# Combine all columns into a single column
combined_data <- waf_data_sub_clean  %>%
  unite("combined", everything(), sep = ",")

# Separate the combined column into individual types
split_data <- combined_data %>%
  separate_rows(combined, sep = ",") %>%
  separate(combined, into = c("Type", "Value"), sep = "_", extra = "merge") %>%
  filter(!is.na(Value))


mut_gene_clean <- split_data %>%
  select(Type, Value) %>%
  group_by(Type, Value) %>%
  summarise(count_gene = n())

# Specify the custom order for Type of Genes
custom_order <- c("E", "M", "N", "NS3", "NS6",  "NS7a",  "NS7b", "NS8", "NSP1", "NSP2", "NSP3", "NSP4", "NSP5", "NSP6", "NSP7", "NSP8", "NSP9", "NSP10", "NSP12", "NSP13", "NSP14", "NSP15", "NSP16", "Spike")
mut_gene_clean$Type <- factor(mut_gene_clean$Type, levels = custom_order)

ggplot(data = mut_gene_clean)+
  geom_bar(aes(x=Type, y= count_gene), stat = "identity", position = "stack")+ 
  labs(x="", y= "Number of mutation") +
  theme_classic() + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 8))
