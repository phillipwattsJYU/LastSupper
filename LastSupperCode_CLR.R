##. PACKAGES

library(phyloseq)
library(qiime2R)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)
library(tibble)
library(scales)
library(stats)
library(forcats)
library(lme4)
library(glmmTMB)
library(nlme)   
library(DHARMa)
library(broom.mixed)
library(performance) 

################################################################################
## DATA PROCESSING
################################################################################
## LOAD QIIME OBJECTS
ps <- qza_to_phyloseq(
  features="TABLE-1-97.qza",
  taxonomy="UNITE-REP-SEQS-1-97.qza",
  metadata = "METADATA-internal.txt" 
)

## access metadata
meta <- sample_data(ps)
str(meta)
# convert to number
meta$timepoint <- as.numeric(as.character(meta$timepoint))
# convert to factor
meta$chip <- factor(as.character(meta$chip))
# assign updated metadata back to the phyloseq object
sample_data(ps) <- meta
# check to ensure the columns are now numeric
str(sample_data(ps))

## CLEAN TAXONOMY TABLE
tax <- data.frame(tax_table(ps))
# remove any leading taxonomy text
tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "d__",""),
                        Phylum = str_replace(tax[,2], "p__",""),
                        Class = str_replace(tax[,3], "c__",""),
                        Order = str_replace(tax[,4], "o__",""),
                        Family = str_replace(tax[,5], "f__",""),
                        Genus = str_replace(tax[,6], "g__",""),
                        Species = str_replace(tax[,7], "s__",""),
                        stringsAsFactors = FALSE)

tax.clean[is.na(tax.clean)] <- ""
# all NAs are now empty cells
for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}

# fill holes in the tax table
tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  # fill in missing taxonomy with last assigned taxonomic level
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}

# insert the cleaned taxonomy table back to phyloseq object
tax_table(ps) <- as.matrix(tax.clean)

# subset to contaminated samples AND remove samples with missing Cs137dose
ps.contam.ad <- subset_samples(ps,
                               treatTime %in% c("Contaminated_1", 
                                                "Contaminated_2") & !is.na(Cs137dose)
)

ps.filt <- prune_samples(sample_sums(ps.contam.ad) >= 400, ps.contam.ad)

## check how much data
ps.filt
microbiome::summarize_phyloseq(ps.filt)

# MACRO / MICRO ration on counts
# load fungal guilds
fungi_class <- read.table(
  "UNITE-REP-SEQS-1-97-fGUILDS.txt",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)
rownames(fungi_class) <- fungi_class$Feature.ID

## match to taxonomy
fungi_class_matched <- fungi_class[rownames(tax_table(ps.filt)), ]

## add Fungal groups to taxonomy
tax.new <- cbind(
  as.data.frame(tax_table(ps.filt)),
  FungiGroup = fungi_class_matched$GrowthMorphology2
)

tax_table(ps.filt) <- tax_table(as.matrix(tax.new))

## aggregate counts by macro / micro categories
ps.mm <- tax_glom(ps.filt, taxrank = "FungiGroup")

## long format
mm_long <- psmelt(ps.mm) %>%
  dplyr::group_by(Sample, FungiGroup) %>%
  dplyr::summarise(Reads = sum(Abundance), .groups = "drop")

## wide format
mm_wide <- mm_long %>%
  tidyr::pivot_wider(
    names_from  = FungiGroup,
    values_from = Reads,
    values_fill = 0
  )

## ensure columns exist
if (!"macro" %in% colnames(mm_wide)) mm_wide$macro <- 0
if (!"micro" %in% colnames(mm_wide)) mm_wide$micro <- 0

## calculate CLR
## add pseudocount
pseudocount <- 0.5

mm_wide <- mm_wide %>%
  dplyr::mutate(
    Ratio_raw = ifelse(micro == 0, Inf, macro / micro),
    CLR_MM = log((macro + pseudocount) / (micro + pseudocount)),
    MM_group = factor(
      ifelse(CLR_MM > 0, "macro", "micro"),
      levels = c("micro", "macro")
    )
  )

## join metadata
meta <- as.data.frame(sample_data(ps.filt))
meta$SampleID <- rownames(meta)

df <- mm_wide %>%
  dplyr::left_join(meta, by = c("Sample" = "SampleID"))

################################################################################
## STATISTICS and PLOTS
################################################################################

# set y limits
max_y <- max(df$Cs137dose, na.rm = TRUE)

## scatterplot
FIG1 <- ggplot(df, aes(x = CLR_MM, y = Cs137dose)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, 
           fill = "#F5FFDC", alpha = 0.3) +  
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, 
           fill = "#8B4513", alpha = 0.3) +  
  geom_point(size = 4, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, 
              linetype = "dashed",  
              color = "black",      
              linewidth = 1.0) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  annotate("text", x = -max(abs(df$CLR_MM), na.rm = TRUE)/2, 
           y = max(df$Cs137dose, na.rm = TRUE)*1.05, 
           label = "more microfungi", hjust = 0.5, size = 6, color = "black") +
  annotate("text", x = max(abs(df$CLR_MM), na.rm = TRUE)/4, 
           y = max(df$Cs137dose, na.rm = TRUE)*1.05, 
           label = "more macrofungi", hjust = 0.5, size = 6, color = "black") +
  labs(
    x = "\nmacrofungi:microfungi ratio (CLR)",
    y = "internal dose of Cs-137 (mGy)\n"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle = 360, hjust = 0.5, size = 18),
    axis.text.y = element_text(angle = 360, hjust = 1, size = 18)
  )

FIG1

## display with one outlier set at limit
FIG1a <- ggplot(df, aes(x = CLR_MM, y = Cs137dose)) +
  annotate("rect", xmin = -6, xmax = 0, ymin = -Inf, ymax = Inf, 
           fill = "#F5FFDC", alpha = 0.3) +  # microfungi-dominated
  annotate("rect", xmin = 0, xmax = 6, ymin = -Inf, ymax = Inf, 
           fill = "#8B4513", alpha = 0.3) +  # macrofungi-dominated
  geom_point(size = 4, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, 
              linetype = "dashed",  
              color = "black",      
              linewidth = 1.0) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  annotate("text", x = -3, 
           y = max(df$Cs137dose, na.rm = TRUE)*1.05, 
           label = "more microfungi", hjust = 0.5, size = 6, color = "black") +
  annotate("text", x = 3, 
           y = max(df$Cs137dose, na.rm = TRUE)*1.05, 
           label = "more macrofungi", hjust = 0.5, size = 6, color = "black") +
  labs(
    x = "\nmacrofungi:microfungi ratio (CLR)",
    y = "internal dose of Cs-137 (mGy)\n"
  ) +
  coord_cartesian(xlim = c(-6, 6)) +  
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle = 360, hjust = 0.5, size = 18),
    axis.text.y = element_text(angle = 360, hjust = 1, size = 18)
  )

FIG1a

## STATISTICAL TEST
## linear model
linear <- lm(
  CLR_MM ~ Cs137dose,
  data = df
)

summary(linear)

## diagnostic plots
# simulate residuals
sim_res <- simulateResiduals(fittedModel = linear, n = 1000)
# plot residual diagnostics
plot(sim_res)

## BOX PLOT
FIG2 <- ggplot(df, aes(x = MM_group, y = Cs137dose, fill = MM_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.8) +   
  geom_jitter(width = 0.15, alpha = 0.7, size = 4) +                 
  scale_fill_manual(values = c("macro" = "#8B4513",
                               "micro" = "#F5FFDC")) +
  labs(
    x = "\ndominant fungal growth form",
    y = "internal dose of Cs-137 (mGy)\n"
  ) +
  scale_y_log10(
    labels = scales::number_format(accuracy = 1),  
    limits = c(0.1, 1000)   
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),                 
    axis.line = element_line(color = "black", linewidth = 0.8),    
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle = 360, hjust = 0.5, size = 18),
    axis.text.y = element_text(angle = 360, hjust = 1, size = 18),
    axis.title.x = element_text(size = 18, face = "plain"),          
    axis.title.y = element_text(size = 18, face = "plain"),         
    legend.position = "none"                                      
  )

FIG2

## Kruskal–Wallis test
kruskal_result <- kruskal.test(Cs137dose ~ MM_group, data = df)
print(kruskal_result)

## calculate median dose per macrofungi category
summary_stats <- df %>%
  group_by(MM_group) %>%
  summarise(
    mean_dose = mean(Cs137dose, na.rm = TRUE),
    median_dose = median(Cs137dose, na.rm = TRUE),
    sd_dose = sd(Cs137dose, na.rm = TRUE),
    n = n()
  )

summary_stats

######################################################
## BARPLOTS of FUNGAL GROWTH FORM by Cs-137 dose
#######################################################
# add traits to tax table
rownames(fungi_class) <- fungi_class$Feature.ID

# extract current taxonomy table as a dataframe 
tax.df <- as.data.frame(tax_table(ps.filt)) %>%
  tibble::rownames_to_column("Feature.ID")

# join fungi_class columns (GrowthMorphology1/2)
tax.extended <- tax.df %>%
  left_join(fungi_class %>% 
              dplyr::select(Feature.ID, GrowthMorphology1, GrowthMorphology2),
            by = "Feature.ID")

# reassign extended taxonomy back to phyloseq
rownames(tax.extended) <- tax.extended$Feature.ID
tax_table(ps.filt) <- tax_table(as.matrix(tax.extended[ , -1]))

# relative abundance
ps.filter.rel <- transform_sample_counts(ps.filt, function(x) x / sum(x))

# agglomerate at GrowthMorphology1 level
glom.gm1 <- tax_glom(ps.filter.rel, taxrank = "GrowthMorphology1", NArm = FALSE)

# convert to long format
data_gm1 <- psmelt(glom.gm1)
data_gm1$GrowthMorphology1 <- as.character(data_gm1$GrowthMorphology1)

# collapse rare groups (<5%)
data_gm1 <- data_gm1 %>%
  dplyr::group_by(Sample, GrowthMorphology1) %>%
  dplyr::summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(GrowthMorphology1 = ifelse(Abundance < 0.05, "<5%", GrowthMorphology1))

# count # phyla to set color palette
Count = length(unique(data_gm1$GrowthMorphology1))
Count # returns 8
unique(data_gm1$GrowthMorphology1)

# order samples by ascending Cs137 dose
meta.df <- data.frame(as(sample_data(ps.filt), "data.frame"))
meta.df$Sample <- rownames(meta.df)

sample_order <- meta.df %>%
  dplyr::arrange(Cs137dose) %>%
  dplyr::pull(Sample)

data_gm1$Sample <- factor(data_gm1$Sample, levels = sample_order)

# define custom order for GrowthMorphology1
gm1_levels <- c("Agaricoid" ,
                "Agaricoid-Boletoid" ,
                "Corticioid" ,
                "Daldinioid" ,
                "Gasteroid" ,
                "Microfungus" ,
                "Yeast",
                "<5%" )   

data_gm1$GrowthMorphology1 <- factor(data_gm1$GrowthMorphology1, levels = gm1_levels)

## make palette
pal <- scales::gradient_n_pal(c("#F5FFDC", "#EAE1C3"))
pal(seq(0, 1, length.out = 4))  # 10 evenly spaced colors

# define custom colors
gm1_colors <- c("Agaricoid" = "#8B4513",
                "Agaricoid-Boletoid" = "#824116",
                "Corticioid" = "#6F3A1A",
                "Daldinioid" = "#8D4123",
                "Gasteroid" = "#66361C",
                "Microfungus" = "#F5FFDC",
                "Yeast" = "#EAE1C3",
                "<5%" = "#CCCCCC")   

# apply same order to both datasets
data_gm1$Sample <- factor(data_gm1$Sample, levels = sample_order)
meta.df$Sample  <- factor(meta.df$Sample,  levels = sample_order)

# growth morphology plot
gm1.plot <- ggplot(data_gm1, aes(x = Sample, y = Abundance, fill = GrowthMorphology1)) +
  geom_bar(stat = "identity", width = 0.9) +
  scale_fill_manual(values = gm1_colors) +
  labs(x = "",
       y = "\nrelative abundance\n",
       fill = "Growth morphology") +   
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),                
    panel.background = element_blank(),          
    axis.line = element_line(color = "black", linewidth = 0.8), 
    axis.ticks.x = element_blank(),    
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.title.x = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    legend.title = element_blank(),        
    legend.text = element_text(size = 14),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

# Cs-137 dose plot
dose.plot <- ggplot(meta.df, aes(x = Sample, y = Cs137dose)) +
  geom_col(fill = "grey50", width = 0.9) +
  scale_y_continuous(
    limits = c(0, 250),                 
    breaks = c(0, 125, 250)   
  ) +
  labs(x = NULL, y = "\nCs-137\n") +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x = element_blank(),   
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "grey50"),
    axis.text.y = element_text(size = 18, color = "grey50"),
    axis.title.y = element_text(size = 18, color = "grey50")
  )

# combine plots
FIG3 <- dose.plot / gm1.plot + plot_layout(heights = c(1, 5))
FIG3

################################################################################
## TEST if MM ratio differs between contaminated and uncontaminated areas
################################################################################

# subset phyloseq object to Contaminated and Uncontaminated
ps.sub <- subset_samples(ps,
                         treatTime %in% c("Contaminated_1", "Contaminated_2",
                                          "Uncontaminated_1", "Uncontaminated_2") &
                           !is.na(Cs137dose))

ps.sub

ps.sub
microbiome::summarize_phyloseq(ps.sub)

## match to taxonomy
fungi_class_matched <- fungi_class[rownames(tax_table(ps.sub)), ]

## add FungiGroup to taxonomy
tax.new <- cbind(
  as.data.frame(tax_table(ps.sub)),
  FungiGroup = fungi_class_matched$GrowthMorphology2
)

tax_table(ps.sub) <- tax_table(as.matrix(tax.new))

## aggregate raw counts by macro / micro
ps.mm <- tax_glom(ps.sub, taxrank = "FungiGroup")

## long format
mm_long <- psmelt(ps.mm) %>%
  dplyr::group_by(Sample, FungiGroup) %>%
  dplyr::summarise(Reads = sum(Abundance), .groups = "drop")

## wide format
mm_wide <- mm_long %>%
  tidyr::pivot_wider(
    names_from  = FungiGroup,
    values_from = Reads,
    values_fill = 0
  )

## ensure columns exist
if (!"macro" %in% colnames(mm_wide)) mm_wide$macro <- 0
if (!"micro" %in% colnames(mm_wide)) mm_wide$micro <- 0

## calculate CLR
pseudocount <- 0.5

mm_wide <- mm_wide %>%
  dplyr::mutate(
    Ratio_raw = ifelse(micro == 0, Inf, macro / micro),
    CLR_MM = log((macro + pseudocount) / (micro + pseudocount)),
    MM_group = factor(
      ifelse(CLR_MM > 0, "macro", "micro"),
      levels = c("micro", "macro")
    )
  )

## join metadata
meta <- as.data.frame(sample_data(ps.sub))
meta$SampleID <- rownames(meta)

df <- mm_wide %>%
  dplyr::left_join(meta, by = c("Sample" = "SampleID"))

# statistical tests
kruskal_result <- kruskal.test(CLR_MM ~ treatment, data = df)
print(kruskal_result)

# boxplot
FIG4 <- ggplot(df, aes(x = treatment, y = CLR_MM, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 4) +
  scale_fill_manual(values = c("Contaminated" = "#8B4513", "Uncontaminated" = "#F5FFDC")) +
  labs(
    x = "\nsample location",
    y = "macrofungi:microfungi ratio (CLR)\n"
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1),
                limits = c(-6, 6), 
                breaks = seq(-6, 6, by = 3)) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.position = "none"
  )

FIG4

###############################################################################################
## Antwis data
###############################################################################################

# LOAD QIIME OBJECTS
ps <- qza_to_phyloseq(
  features="TABLE-97-EE8-group.qza",
  taxonomy="UNITE-REP-SEQS-97-EE8.qza",
  metadata = "METADATA-group-1.txt"
)

## CLEAN TAXONOMY TABLE
tax <- data.frame(tax_table(ps))
tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "d__",""),
                        Phylum = str_replace(tax[,2], "p__",""),
                        Class = str_replace(tax[,3], "c__",""),
                        Order = str_replace(tax[,4], "o__",""),
                        Family = str_replace(tax[,5], "f__",""),
                        Genus = str_replace(tax[,6], "g__",""),
                        Species = str_replace(tax[,7], "s__",""),
                        stringsAsFactors = FALSE)

tax.clean[is.na(tax.clean)] <- ""
for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}

tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}

tax_table(ps) <- as.matrix(tax.clean)

# subset to data - remove NAs and keep Clethrionomys
ps.data <- subset_samples(ps, species == "Myodes_glareolus")

# subset to contaminated area
ps.contam <- subset_samples(ps.data, treatment == "cont")

ps.filt <- prune_samples(sample_sums(ps.contam) >= 100, ps.contam)

fungi_class <- read.table(
  "UNITE-REP-SEQS-97-EE8-fGUILD.txt",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

rownames(fungi_class) <- fungi_class$Feature.ID

## match to taxonomy
fungi_class_matched <- fungi_class[rownames(tax_table(ps.filt)), ]

## add FungiGroup to taxonomy
tax.new <- cbind(
  as.data.frame(tax_table(ps.filt)),
  FungiGroup = fungi_class_matched$GrowthMorphology2
)

tax_table(ps.filt) <- tax_table(as.matrix(tax.new))

## aggregate raw counts by macro / micro
ps.mm <- tax_glom(ps.filt, taxrank = "FungiGroup")

## long format
mm_long <- psmelt(ps.mm) %>%
  dplyr::group_by(Sample, FungiGroup) %>%
  dplyr::summarise(Reads = sum(Abundance), .groups = "drop")

## wide format
mm_wide <- mm_long %>%
  tidyr::pivot_wider(
    names_from  = FungiGroup,
    values_from = Reads,
    values_fill = 0
  )

## ensure columns exist
if (!"macro" %in% colnames(mm_wide)) mm_wide$macro <- 0
if (!"micro" %in% colnames(mm_wide)) mm_wide$micro <- 0

## CLR
pseudocount <- 0.5

mm_wide <- mm_wide %>%
  dplyr::mutate(
    # protect against division by zero
    Ratio_raw = ifelse(micro == 0, Inf, macro / micro),
    CLR_MM = log((macro + pseudocount) / (micro + pseudocount)),
    # categorical groupings based on CLR
    MM_group = factor(
      ifelse(CLR_MM > 0, "macro", "micro"),
      levels = c("micro", "macro")
    )
  )

## join metadata
meta <- as.data.frame(sample_data(ps.filt))
meta$SampleID <- rownames(meta)

df <- mm_wide %>%
  dplyr::left_join(meta, by = c("Sample" = "SampleID"))

## check data
ps.filt
microbiome::summarize_phyloseq(ps.filt)

################################################################################
## STATISTICS and PLOTS
################################################################################

## scatterplot
FIG1 <- ggplot(df, aes(x = CLR_MM, y = cs137dose)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, 
           fill = "#F5FFDC", alpha = 0.3) +  
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, 
           fill = "#8B4513", alpha = 0.3) +  
  geom_point(size = 4, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, 
              linetype = "dashed",  
              color = "black",      
              linewidth = 1.0) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  annotate("text", x = -max(abs(df$CLR_MM), na.rm = TRUE)/2, 
           y = max(df$cs137dose, na.rm = TRUE)*1.05, 
           label = "more microfungi", hjust = 0.5, size = 6, color = "black") +
  annotate("text", x = max(abs(df$CLR_MM), na.rm = TRUE)/2, 
           y = max(df$cs137dose, na.rm = TRUE)*1.05, 
           label = "more macrofungi", hjust = 0.5, size = 6, color = "black") +
  labs(
    x = "\nmacrofungi:microfungi ratio (CLR)",
    y = "internal dose of Cs-137 (mGy)\n"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle = 360, hjust = 0.5, size = 18),
    axis.text.y = element_text(angle = 360, hjust = 1, size = 18)
  )

FIG1

## STATISTICAL TEST - MM versus CLR
linear <- lm(
  CLR_MM ~ cs137dose,
  data = df
)

summary(linear)

# simulate residuals
sim_res <- simulateResiduals(fittedModel = linear, n = 1000)
# plot residual diagnostics
plot(sim_res)

## test groups
## Kruskal–Wallis test
kruskal_result <- kruskal.test(cs137dose ~ MM_group, data = df)
print(kruskal_result)

## calculate median dose per category
summary_stats <- df %>%
  group_by(MM_group) %>%
  summarise(
    mean_dose = mean(cs137dose, na.rm = TRUE),
    median_dose = median(cs137dose, na.rm = TRUE),
    sd_dose = sd(cs137dose, na.rm = TRUE),
    n = n()
  )

summary_stats

## BOXPLOTS
FIG2 <- ggplot(df, aes(x = MM_group, y = cs137dose, fill = MM_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.8) +   
  geom_jitter(width = 0.15, alpha = 0.7, size = 4) +                 
  scale_fill_manual(values = c("macro" = "#8B4513",
                               "micro" = "#F5FFDC")) +
  labs(
    x = "\ndominant fungal growth form",
    y = "internal dose of Cs-137 (mGy)\n"
  ) +
  scale_y_log10(
    labels = scales::number_format(accuracy = 1),  
    limits = c(0.1, 1000)   
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),                 
    axis.line = element_line(color = "black", linewidth = 0.8),    
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle = 360, hjust = 0.5, size = 18),
    axis.text.y = element_text(angle = 360, hjust = 1, size = 18),
    axis.title.x = element_text(size = 18, face = "plain"),          
    axis.title.y = element_text(size = 18, face = "plain"),         
    legend.position = "none"                                      
  )

FIG2

######################################################
## BARPLOTS of growth form
#######################################################

# add traits to tax table
rownames(fungi_class) <- fungi_class$Feature.ID

# extract current taxonomy table as a dataframe 
tax.df <- as.data.frame(tax_table(ps.filt)) %>%
  tibble::rownames_to_column("Feature.ID")

# join fungi_class columns (GrowthMorphology1/2)
tax.extended <- tax.df %>%
  left_join(fungi_class %>% 
              dplyr::select(Feature.ID, GrowthMorphology1, GrowthMorphology2),
            by = "Feature.ID")

# reassign extended taxonomy back to phyloseq
rownames(tax.extended) <- tax.extended$Feature.ID
tax_table(ps.filt) <- tax_table(as.matrix(tax.extended[ , -1]))

# relative abundance
ps.filter.rel <- transform_sample_counts(ps.filt, function(x) x / sum(x))

# agglomerate at GrowthMorphology1 level
glom.gm1 <- tax_glom(ps.filter.rel, taxrank = "GrowthMorphology1", NArm = FALSE)

# convert to long format
data_gm1 <- psmelt(glom.gm1)
data_gm1$GrowthMorphology1 <- as.character(data_gm1$GrowthMorphology1)

# collapse rare groups (<5%)
data_gm1 <- data_gm1 %>%
  dplyr::group_by(Sample, GrowthMorphology1) %>%
  dplyr::summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(GrowthMorphology1 = ifelse(Abundance < 0.05, "<5%", GrowthMorphology1))

# count # phyla to set color palette
Count = length(unique(data_gm1$GrowthMorphology1))
Count # returns 11
unique(data_gm1$GrowthMorphology1)

# order samples by ascending Cs137 dose
meta.df <- data.frame(as(sample_data(ps.filt), "data.frame"))
meta.df$Sample <- rownames(meta.df)

sample_order <- meta.df %>%
  dplyr::arrange(cs137dose) %>%
  dplyr::pull(Sample)

data_gm1$Sample <- factor(data_gm1$Sample, levels = sample_order)

# define custom order for GrowthMorphology1
gm1_levels <- c("Agaricoid",
                "Agaricoid-Boletoid",
                "Clavarioid",
                "Corticioid",
                "Gasteroid",
                "Pezizoid",
                "Polyporoid",
                "Macrofungus",
                "Microfungus",
                "Yeast",
                "<5%")   

data_gm1$GrowthMorphology1 <- factor(data_gm1$GrowthMorphology1, levels = gm1_levels)

## make palette
pal <- scales::gradient_n_pal(c("#8B4513", "#4B2B1F"))
pal(seq(0, 1, length.out = 8))  

# define custom colors
gm1_colors <- c("Agaricoid" = "#CF4616",
                "Agaricoid-Boletoid" = "#BE451A",
                "Clavarioid" = "#AD441E",
                "Corticioid" = "#9D4321",
                "Gasteroid" = "#8D4123",
                "Pezizoid" = "#7C3F25",
                "Polyporoid" = "#6C3D27",
                "Macrofungus" = "#5C3A29",
                "Microfungus" = "#F5FFDC",
                "Yeast" = "#EAE1C3",
                "<5%" = "#CCCCCC")  

gm1_colors <- c("Agaricoid" = "#8B4513",
                "Agaricoid-Boletoid" = "#824116",
                "Clavarioid" = "#783D18",
                "Corticioid" = "#6F3A1A",
                "Gasteroid" = "#8D4123",
                "Pezizoid" = "#66361C",
                "Polyporoid" = "#542F1E",
                "Macrofungus" = "#4B2B1F",
                "Microfungus" = "#F5FFDC",
                "Yeast" = "#EAE1C3",
                "<5%" = "#CCCCCC")  

# apply same order to both datasets
data_gm1$Sample <- factor(data_gm1$Sample, levels = sample_order)
meta.df$Sample  <- factor(meta.df$Sample,  levels = sample_order)

# growth morphology plot
gm1.plot <- ggplot(data_gm1, aes(x = Sample, y = Abundance, fill = GrowthMorphology1)) +
  geom_bar(stat = "identity", width = 0.9) +
  scale_fill_manual(values = gm1_colors) +
  labs(x = "",
       y = "\nrelative abundance\n",
       fill = "Growth morphology") +  
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),             
    panel.background = element_blank(),     
    axis.line = element_line(color = "black", linewidth = 0.8), # black axes
    axis.ticks.x = element_blank(),    
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.title.x = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    legend.title = element_blank(),      
    legend.text = element_text(size = 14),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

# Cs-137 dose plot
dose.plot <- ggplot(meta.df, aes(x = Sample, y = cs137dose)) +
  geom_col(fill = "grey50", width = 0.9) +
  #  scale_y_log10(
  scale_y_continuous(
    limits = c(0, 1000),      
    breaks = c(0, 500, 1000)   
  ) +
  labs(x = NULL, y = "\nCs-137\n") +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x = element_blank(),   
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "grey50"),
    axis.text.y = element_text(size = 18, color = "grey50"),
    axis.title.y = element_text(size = 18, color = "grey50")
  )

# combine plots
FIG3 <- dose.plot / gm1.plot + plot_layout(heights = c(1, 5))
FIG3

################################################################################
## TEST if MM ratio differs between cont and uncont areas
################################################################################
# subset to Contaminated and Uncontaminated
ps.sub <- subset_samples(ps,
                         treatment %in% c("cont", "uncon") &
                           species1 %in% c("C.gla") &
                           !is.na(cs137dose))

ps.sub.filter <- prune_samples(
  sample_sums(ps.sub) > 100 & sample_data(ps.sub)$species1 == "C.gla",
  ps.sub
)

ps.sub.filter
microbiome::summarize_phyloseq(ps.sub.filter)

fungi_class <- read.table(
  "UNITE-REP-SEQS-97-EE8-fGUILD.txt",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

rownames(fungi_class) <- fungi_class$Feature.ID

## match to taxonomy
fungi_class_matched <- fungi_class[rownames(tax_table(ps.sub.filter)), ]

## add FungiGroup to taxonomy
tax.new <- cbind(
  as.data.frame(tax_table(ps.sub.filter)),
  FungiGroup = fungi_class_matched$GrowthMorphology2
)

tax_table(ps.sub.filter) <- tax_table(as.matrix(tax.new))

## aggregate RAW counts by macro / micro
ps.mm <- tax_glom(ps.sub.filter, taxrank = "FungiGroup")

## long format (Sample × Group × Reads)
mm_long <- psmelt(ps.mm) %>%
  dplyr::group_by(Sample, FungiGroup) %>%
  dplyr::summarise(Reads = sum(Abundance), .groups = "drop")

## wide format
mm_wide <- mm_long %>%
  tidyr::pivot_wider(
    names_from  = FungiGroup,
    values_from = Reads,
    values_fill = 0
  )

## ensure columns exist
if (!"macro" %in% colnames(mm_wide)) mm_wide$macro <- 0
if (!"micro" %in% colnames(mm_wide)) mm_wide$micro <- 0

## CLR
pseudocount <- 0.5

mm_wide <- mm_wide %>%
  dplyr::mutate(
    Ratio_raw = ifelse(micro == 0, Inf, macro / micro),
    CLR_MM = log((macro + pseudocount) / (micro + pseudocount)),
    MM_group = factor(
      ifelse(CLR_MM > 0, "macro", "micro"),
      levels = c("micro", "macro")
    )
  )

## join metadata
meta <- as.data.frame(sample_data(ps.sub.filter))
meta$SampleID <- rownames(meta)

df <- mm_wide %>%
  dplyr::left_join(meta, by = c("Sample" = "SampleID"))

# statistical tests
kruskal_result <- kruskal.test(CLR_MM ~ treatment, data = df)
print(kruskal_result)

# boxplot (contaminated vs uncontaminated pooled)
FIG4 <- ggplot(df, aes(x = treatment, y = CLR_MM, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 4) +
  scale_fill_manual(values = c("cont" = "#8B4513", "uncon" = "#F5FFDC")) +
  labs(
    x = "\nsample location",
    y = "macrofungi:microfungi ratio (CLR)\n"
  ) +
  scale_x_discrete(
    labels = c(
      "cont"  = "Contaminated",
      "uncon" = "Uncontaminated"
    )
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1),
                limits = c(-12, 12), 
                breaks = seq(-12, 12, by = 6)) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.position = "none"
  )

FIG4
