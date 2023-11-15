# Group Data Analysis Project

This group assignment was a whole-term data analysis project analysing cancer genomics data obtained from the [COSMIC](https://cancer.sanger.ac.uk/cosmic) database. 

The datasets that we analysed were: cancer samples, cancer hallmarks, cancer gene expression, cancer mutations and cancer resistance mutations.

Different group members were more (or less) heavily involved in the analysis of different datasets, and in my case, I was primarily invovled in analysing gene expression, mutations and resistance mutations datasets. I will highlight my personal contributions to the project below, but the whole, finished project with everyone's contributions and notes can be accessed in the *Buck Compiled Code.Rmd* file. 

## Cancer gene expression and mutations
I started by exploring the cancer gene expression datased and performing the regulation count, plotting the scatterplot demonstrating that over-expressed genes with highest z-scores are outliers in the data, and coding for the mean and SD of over-, under- and normally regulated genes (including the outliers).  
```
a <- gene_expression %>%
  count(Regulation, sort = TRUE) %>%
  mutate(props = prop.table(n))
a

gene_expression %>%
  select(`Gene Name`, Regulation, `Z Score`) %>%
  filter(Regulation == "over") %>%
  ggplot() + aes(x = `Gene Name`, y = `Z Score`) + geom_point() + theme(axis.text.x = element_text(angle=90, hjust=1, size = 1)) + labs(title = "Over Expressed genes and their respective Z scores")

gene_expression %>%
  select(`Gene Name`, `Regulation`, `Z Score`) %>%
  filter(Regulation == "over", `Z Score` > 58) %>%
  ggplot() + aes(x = `Gene Name`, y = `Z Score`) + geom_point() + theme(axis.text.x = element_text(angle=90, hjust=1)) + labs(title = "Genes with incredibly high Z scores") + labs(subtitle = "(Z Score > 58)")

# Mean and standard deviations of z-scores for over-expressed genes
gene_expression %>%
  filter(`Z Score` > 2) %>%
  arrange(`Z Score`) %>%
    summarise(across(.cols = c(`Z Score`), list(average = ~mean(.x, na.rm = TRUE), stdev = ~sd(.x, na.rm = TRUE))))

# Mean and standard deviations of z-scores for normally regulated genes
gene_expression %>%
  select(Regulation, `Z Score`) %>%
  filter(Regulation == "normal") %>%
  summarise(across(.cols = c(`Z Score`), list(average = ~mean(.x, na.rm = TRUE), stdev = ~sd(.x, na.rm = TRUE))))

# Mean and standard deviations of z-scores for under-expressed genes
gene_expression %>%
  select(Regulation, `Z Score`) %>%
  filter(Regulation == "under") %>%
  summarise(across(.cols = c(`Z Score`), list(average = ~mean(.x, na.rm = TRUE), stdev = ~sd(.x, na.rm = TRUE))))
```

I have then merged the cancer gene expression and cancer mutations datasets by gene name to create over-expressed and under-expressed data frames, performed a count of over- and under- expressed genes, primary histologies, primary histologies of the most frequent over- and under- expressed gene, and genes over- and under- expressed in the most frequent primary histologies in both of the new data frames. I've included the code excerpt from this exploratory analysis of the over-expressed data frame below (the same code was used to analyse the under-expressed data frame). 
```
# Find top over-expressed genes in cancer
overexpressed <- inner_join(gene_expression, mutation, by = "Gene Name") %>%
  select(`Gene Name`, Regulation, `Z Score`, `Primary Site`, `Primary Histology`, `Histology Subtype 1`, `Fathmm Prediction`, `Fathmm Score`, Age,) %>%
  filter(Regulation == "over", `Fathmm Prediction` == "PATHOGENIC") %>%
  arrange(desc(`Fathmm Score`))
overexpressed

overexpressed %>%
  count(`Gene Name`) %>%
  arrange(desc(n))

# Count of primary histologies
overexpressed %>%
  count(`Primary Histology`) %>%
  arrange(desc(n))

overexpressed %>%
  count(`Primary Histology`) %>%
  filter(n > 8) %>%
  arrange(desc(n))

# Count of genes over-expressed in the most frequent primary histology
overexpressed %>%
  filter(`Primary Histology` == "carcinoma") %>%
  count(`Gene Name`) %>%
  filter(n > 7) %>%
  arrange(desc(n)) 
```

My colleague and I have then filtered the top over- and under-expressed genes, plotted a graph demonstrating a high confidence of MTOR being under- and TP53 being over-regulated in general and in carcinoma (most frequent primary histology), and plotted a boxplot of 4 top primary histologies against fathmm scores of top over- and under-expressed genes. 
```
# Primary histologies of the most frequently over-expressed gene (TP53)
overexpressed %>%
  filter(`Gene Name` == "TP53") %>%
  count(`Primary Histology`) %>%
  arrange(desc(n)) %>%
  mutate(`Primary Histology` = case_when(
    `Primary Histology` == "adnexal_tumour" ~ "Adnexal Tumour",
    `Primary Histology` == "adnexal_cortical_neoplasm" ~ "Adnexal Cortical Neoplasm",
    `Primary Histology` == "carcinoma" ~ "Carcinoma",
    `Primary Histology` == "germ_cell_tumour" ~ "Germ Cell Tumour",
    `Primary Histology` == "glioma" ~ "Glioma",
    `Primary Histology` == "haematopoietic_neoplasm" ~ "Haematopoietic Neoplasm",
    `Primary Histology` == "leiomyosarcoma" ~ "Leiomyosarcoma",
    `Primary Histology` == "lymphoid_neoplasm" ~ "Lymphoid Neoplasm",
    `Primary Histology` == "malignant_melanoma" ~ "Malignant Melanoma",
    `Primary Histology` == "osteosarcoma" ~ "Osteosarcoma",
    `Primary Histology` == "other" ~ "Other", 
    `Primary Histology` == "thymoma" ~ "Thymoma",
    `Primary Histology` == "Wilms_tumour" ~ "Wilms Tumour")) %>%
  drop_na(`Primary Histology`) %>%
  ggplot() + aes(x = `Primary Histology`, y = n, fill = `Primary Histology`) + geom_col() + ylab("Frequency") + labs(title = " Frequency of different Primary histologies") + labs(subtitle = "TP53 Gene") + theme(axis.text.x = element_text(angle = 90))

FATHMM_Boxplot <- inner_join(gene_expression, mutation, by = "Gene Name") %>%
  select(`Gene Name`, `Primary Histology`, `Fathmm Prediction`, `Fathmm Score`) %>%
  filter(`Gene Name` %in% c("MTOR", "KDM5C", "RB1", "NSD1", "TP53", "EGFR", "KMT2C", "FAT3", "NOTCH1", "FBXW7", "ATR", "LRIG3", "GRIN2A"), `Fathmm Prediction` == "PATHOGENIC", `Primary Histology` %in% c("carcinoma", "malignant_melanoma", "glioma", "lymphoid_neoplasm")) %>%
  mutate(`Primary Histology` = case_when(
    `Primary Histology`== "carcinoma" ~ "Carcinoma",
    `Primary Histology`== "glioma" ~ "Glioma",
    `Primary Histology`== "lymphoid_neoplasm" ~ "Lymphoid Neoplasm",
    `Primary Histology`== "malignant_melanoma" ~ "Malignant Melanoma")) %>%
  ggplot() + aes(x = `Primary Histology`, y = `Fathmm Score`, colour = `Gene Name`) + geom_boxplot() +
  labs(title = "Boxplot of FATHMM Scores of most commonly over- and under-expressed genes") +
  theme(plot.title = element_text(size = rel(1))) +
  xlab("Primary Histology") + ylab("FATHMM Score") + labs(subtitle = "Top 4 most common Primary Histologies")
```

## Cancer gene expression and resistance mutations
After finishing the first part of the project which mainly involved exploring gene expression and mutations, we now moved on to look into the pharmacological side of things and drug resistance mutations in cancer. First, I merged the gene expression and resistance mutations datasets, and then I performed gene name, regulation, histology, drug name and tier counts.
```
# Drug resistance analysis
drugs <- inner_join(gene_expression, r_mutation, by = "Gene Name")
drugs %>%
  count(`Gene Name`) %>%
  arrange(desc(n))

drugs %>%
  count(`Regulation`)

drugs %>%
  count(`Histology`) %>%
  arrange(desc(n))

drugs %>%
  count(`Drug Name`) %>%
  arrange(desc(n))

drugs %>%
  count(`Tier`)
```

I have plotted a drug name vs histology graph, used `facet_wrap` function to divide it by gene names, and coloured it by primary tissue, demonstrating which genes show resistance to which drugs, in which histologies and in which primary tissues.  
```
drugs %>%
  select(`Gene Name`, `Drug Name`, `Regulation`, `Primary Tissue`, `Histology`, `Histology Subtype 1`) %>%
  filter(`Regulation` != "normal", `Histology` != c("epithelioid_inflammatory_myofibroblastic_sarcoma", "inflammatory_myofibroblastic_tumour", "primitive_neuroectodermal_tumour-medulloblastoma")) %>%
  distinct() %>%
  mutate(`Primary Tissue` = case_when(
    `Primary Tissue`== "breast" ~ "Breast",
    `Primary Tissue`== "haematopoietic_and_lymphoid_tissue" ~ "Haematopoietic and Lymphoid Tissue",
    `Primary Tissue`== "kidney" ~ "Kidney",
    `Primary Tissue`== "lung" ~ "Lung",
    `Primary Tissue` == "thyroid" ~ "Thyroid")) %>%
    mutate(`Histology` = case_when(
    `Histology` == "carcinoma" ~ "Carcinoma",
    `Histology` == "haematopoietic_neoplasm" ~ "Haematopoietic Neoplasm",
    `Histology` == "lymphoid_neoplasm" ~ "Lymphoid Neoplasm")) %>%
  ggplot() + aes(x = `Histology`, y = `Drug Name`, colour = `Primary Tissue`) + geom_point() + theme(axis.text.x = element_text(size = 8.5, angle = 90)) + facet_wrap(~`Gene Name`)
```

Last part of the project involved incorporating information from the hallmarks dataset into the gene expression and resistance mutation analysis, so I merged hallmarks and drugs data frames, and mapped hallmarks, genome coordinates and impact to EGFR gene in lung carcinoma.  
```
celltypes <- inner_join(drugs, hallmarks, by = "Gene Name")

celltypes %>%
  select(`Gene Name`, `Drug Name`, `Cell Type`, `Hallmark`, `Impact`, `Regulation`, `Primary Tissue`, `Histology`, `Histology Subtype 1`, `Pubmed Id`) %>%
  filter(`Regulation` != "normal", `Primary Tissue` == "lung", `Histology Subtype 1` %in% c("adenocarcinoma", "non_small_cell_carcinoma"), `Gene Name` == "EGFR") %>% 
  count(Hallmark) %>%
  arrange(desc(n))

celltypes %>%
  select(`Gene Name`, `Drug Name`, `Cell Type`, `Hallmark`, `Impact`, `Regulation`, `Primary Tissue`, `Histology`, `Histology Subtype 1`, `Genome Coordinates Gr Ch38`, `Pubmed Id`) %>%
  filter(`Regulation` != "normal", `Primary Tissue` == "lung", `Histology Subtype 1` %in% c("adenocarcinoma", "non_small_cell_carcinoma"), `Gene Name` == "EGFR") %>% 
  drop_na(`Genome Coordinates Gr Ch38`) %>%
  count(Hallmark, `Genome Coordinates Gr Ch38`, Impact, `Pubmed Id`) %>%
  arrange(desc(n))
```
