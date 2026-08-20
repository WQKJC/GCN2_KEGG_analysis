
##-----------Virtual Knockout------------------##

library(dplyr)
library(ggplot2)



my_selected_terms <- c(
  "Cytoskeleton in muscle cells",
  "Focal adhesion",
  "Regulation of actin cytoskeleton",
  "ECM-receptor interaction",
  "Vascular smooth muscle contraction",
  "Adherens junction",
  "PI3K-Akt signaling pathway",
  "MAPK signaling pathway",
  "TNF signaling pathway",
  "mTOR signaling pathway"
)


kegg_data <- read.csv("enrichment-kegg-All-vs-DEG-Total-Virtual Knockout .csv", header = TRUE, stringsAsFactors = FALSE, encoding = "UTF-8") %>% 
  filter(Description %in% my_selected_terms)

kegg_data$Description <- factor(kegg_data$Description, levels = rev(my_selected_terms))


p <- ggplot(kegg_data, aes(x = NES, y = Description)) +
  geom_point(aes(size = setSize, color = pvalue)) +
  

  scale_color_gradientn(
    colors = c("#b2182b", "#fdae61", "#f7f7f7", "#abd9e9", "#2166ac"),
    limits = c(0, 0.05) 
  ) + 
 
  scale_size_continuous(
    range = c(4, 11),
    limits = c(0, 200)  
  ) +
  
  theme_bw() + 
  labs(
    title = "KEGG analysis",
    x = "Enrichment Score",
    y = NULL,
    color = "p_value",
    size = "Gene Count"
  ) +
  
  guides(
    color = guide_colorbar(
      title.position = "top",
      barwidth = unit(1.2, "lines"), 
      barheight = unit(6, "lines"),  
      order = 2
    ),
    size = guide_legend(
      title.position = "top",
      override.aes = list(color = "black"),
      order = 1
    )
  ) +
  

  theme(
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"), 
    legend.position = "right",             
    legend.box = "vertical",               
    legend.margin = margin(l = 15, unit = "pt"), 
    legend.spacing.y = unit(20, "pt"),           
    legend.title = element_text(face = "bold", size = 11, color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(face = "bold", color = "black", size = 13), 
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b=15)),
    
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray95")
  )

pdf("Virtual Knockout.pdf", width = 8.5, height = 8)
print(p)
dev.off()




##-----------lung-MMC------------------##

library(dplyr)
library(ggplot2)



my_selected_terms <- c(
  "ECM-receptor interaction",
  "Complement and coagulation cascades",
  "Focal adhesion",
  "ErbB signaling pathway",
  "PI3K-Akt signaling pathway",
  "Diabetic cardiomyopathy",
  "Fluid shear stress and atherosclerosis",
  "Pathways in cancer",
  "MAPK signaling pathway",
  "Cytokine-cytokine receptor interaction"
)


kegg_raw <- read.csv("enrichment-kegg-All-vs-DEG-Total-lung-MMC.csv", header = TRUE, stringsAsFactors = FALSE, encoding = "UTF-8") %>% 
  filter(Term %in% my_selected_terms)

colnames(kegg_raw)[colnames(kegg_raw) == "p.value"] <- "p_value"

kegg_data <- kegg_raw %>% 
  filter(Term %in% my_selected_terms)
kegg_data$Term <- factor(kegg_data$Term, levels = rev(my_selected_terms))

p <- ggplot(kegg_data, aes(x = Enrichment_score, y = Term)) +
  geom_point(aes(size = ListHits, color = p_value)) +
  
  
  scale_color_gradientn(
    colors = c("#b2182b", "#fdae61", "#f7f7f7", "#abd9e9", "#2166ac"),
    limits = c(0, 0.05) 
  ) + 
  
  scale_size_continuous(
    range = c(4, 11),
    limits = c(0, 200)  
  ) +
  
  theme_bw() + 
  labs(
    title = "KEGG analysis",
    x = "Enrichment Score",
    y = NULL,
    color = "p_value",
    size = "Gene Count"
  ) +
  
  guides(
    color = guide_colorbar(
      title.position = "top",
      barwidth = unit(1.2, "lines"), 
      barheight = unit(6, "lines"),  
      order = 2
    ),
    size = guide_legend(
      title.position = "top",
      override.aes = list(color = "black"),
      order = 1
    )
  ) +
  
  
  theme(
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"), 
    legend.position = "right",             
    legend.box = "vertical",               
    legend.margin = margin(l = 15, unit = "pt"), 
    legend.spacing.y = unit(20, "pt"),           
    legend.title = element_text(face = "bold", size = 11, color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(face = "bold", color = "black", size = 13), 
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b=15)),
    
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray95")
  )

pdf("lung-MMC.pdf", width = 8.5, height = 8)
print(p)
dev.off()



##-----------lung-WT------------------##

library(dplyr)
library(ggplot2)



my_selected_terms <- c(
  "TNF signaling pathway",
  "NF-kappa B signaling pathway",
  "Acute myeloid leukemia",
  "Cytokine-cytokine receptor interaction",
  "IL-17 signaling pathway",
  "Fluid shear stress and atherosclerosis",
  "Cell adhesion molecules",
  "Chemokine signaling pathway",
  "Phagosome",
  "MAPK signaling pathway"
  
)


kegg_raw <- read.csv("enrichment-kegg-All-vs-DEG-Total-lung.csv", header = TRUE, stringsAsFactors = FALSE, encoding = "UTF-8") %>% 
  filter(Term %in% my_selected_terms)



colnames(kegg_raw)[colnames(kegg_raw) == "p.value"] <- "p_value"


kegg_data <- kegg_raw %>% 
  filter(Term %in% my_selected_terms)

kegg_data$Term <- factor(kegg_data$Term, levels = rev(my_selected_terms))


p <- ggplot(kegg_data, aes(x = Enrichment_score, y = Term)) +
  geom_point(aes(size = ListHits, color = p_value)) +
  

  scale_color_gradientn(
    colors = c("#b2182b", "#fdae61", "#f7f7f7", "#abd9e9", "#2166ac"),
    limits = c(0, 0.05) 
  ) + 
  
  
  scale_size_continuous(
    range = c(5, 11),
    limits = c(0, 200)  
  ) +
  
  theme_bw() + 
  labs(
    title = "KEGG analysis",
    x = "Enrichment Score",
    y = NULL,
    color = "p_value",
    size = "Gene Count"
  ) +
  
 
  guides(
    color = guide_colorbar(
      title.position = "top",
      barwidth = unit(1.2, "lines"), 
      barheight = unit(6, "lines"),  
      order = 2
    ),
    size = guide_legend(
      title.position = "top",
      override.aes = list(color = "black"),
      order = 1
    )
  ) +
  
  
  theme(
    
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"), 
    
    
    legend.position = "right",             
    legend.box = "vertical",               
    legend.margin = margin(l = 15, unit = "pt"), 
    legend.spacing.y = unit(20, "pt"),           
    
    
    legend.title = element_text(face = "bold", size = 11, color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(face = "bold", color = "black", size = 13), 
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b=15)),
    
    
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray95")
  )


pdf("lung-WT.pdf", width = 8.5, height = 8)
print(p)
dev.off()




##-----------MMC-SMC------------------##



library(dplyr)
library(ggplot2)



my_selected_terms <- c(
  "Long-term potentiation",
  "Longevity regulating pathway",
  "Hematopoietic cell lineage",
  "Chemokine signaling pathway",
  "cAMP signaling pathway",
  "Cholinergic synapse",
  "Oocyte meiosis",
  "Apelin signaling pathway",
  "Vascular smooth muscle contraction",
  "Cytokine-cytokine receptor interaction"
)


kegg_raw <- read.csv("enrichment-kegg-All-vs-DEG-Total-MMC-SMC.csv", header = TRUE, stringsAsFactors = FALSE, encoding = "UTF-8") %>% 
  filter(Term %in% my_selected_terms)



colnames(kegg_raw)[colnames(kegg_raw) == "p.value"] <- "p_value"


kegg_data <- kegg_raw %>% 
  filter(Term %in% my_selected_terms)

kegg_data$Term <- factor(kegg_data$Term, levels = rev(my_selected_terms))


p <- ggplot(kegg_data, aes(x = Enrichment_score, y = Term)) +
  geom_point(aes(size = ListHits, color = p_value)) +
  
  
  scale_color_gradientn(
    colors = c("#b2182b", "#fdae61", "#f7f7f7", "#abd9e9", "#2166ac"),
    limits = c(0, 0.05) 
  ) + 
  
  
  scale_size_continuous(
    range = c(5, 11),
    limits = c(0, 200)  
  ) +
  
  theme_bw() + 
  labs(
    title = "KEGG analysis",
    x = "Enrichment Score",
    y = NULL,
    color = "p_value",
    size = "Gene Count"
  ) +
  
  
  guides(
    color = guide_colorbar(
      title.position = "top",
      barwidth = unit(1.2, "lines"), 
      barheight = unit(6, "lines"),  
      order = 2
    ),
    size = guide_legend(
      title.position = "top",
      override.aes = list(color = "black"),
      order = 1
    )
  ) +
  
  
  theme(
    
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"), 
    
    
    legend.position = "right",             
    legend.box = "vertical",               
    legend.margin = margin(l = 15, unit = "pt"), 
    legend.spacing.y = unit(20, "pt"),           
    
    
    legend.title = element_text(face = "bold", size = 11, color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(face = "bold", color = "black", size = 13), 
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b=15)),
    
    
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray95")
  )

pdf("MMC-SMC.pdf", width = 8.5, height = 8)
print(p)
dev.off()




##-----------WT-SMC------------------##


library(dplyr)
library(ggplot2)



my_selected_terms <- c(
  "ECM-receptor interaction",
  "IL-17 signaling pathway",
  "Fluid shear stress and atherosclerosis",
  "Cell adhesion molecules",
  "B cell receptor signaling pathway",
  "TNF signaling pathway",
  "Leukocyte transendothelial migration",
  "PI3K-Akt signaling pathway",
  "Cellular senescence",
  "MAPK signaling pathway"
)


kegg_raw <- read.csv("enrichment-kegg-All-vs-DEG-Total-SMC.csv", header = TRUE, stringsAsFactors = FALSE, encoding = "UTF-8") %>% 
  filter(Term %in% my_selected_terms)



colnames(kegg_raw)[colnames(kegg_raw) == "p.value"] <- "p_value"


kegg_data <- kegg_raw %>% 
  filter(Term %in% my_selected_terms)

kegg_data$Term <- factor(kegg_data$Term, levels = rev(my_selected_terms))


p <- ggplot(kegg_data, aes(x = Enrichment_score, y = Term)) +
  geom_point(aes(size = ListHits, color = p_value)) +
  
  
  scale_color_gradientn(
    colors = c("#b2182b", "#fdae61", "#f7f7f7", "#abd9e9", "#2166ac"),
    limits = c(0, 0.05) 
  ) + 
  
  
  scale_size_continuous(
    range = c(5, 11),
    limits = c(0, 200)  
  ) +
  
  theme_bw() + 
  labs(
    title = "KEGG analysis",
    x = "Enrichment Score",
    y = NULL,
    color = "p_value",
    size = "Gene Count"
  ) +
  
  
  guides(
    color = guide_colorbar(
      title.position = "top",
      barwidth = unit(1.2, "lines"), 
      barheight = unit(6, "lines"),  
      order = 2
    ),
    size = guide_legend(
      title.position = "top",
      override.aes = list(color = "black"),
      order = 1
    )
  ) +
  
  
  theme(
   
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"), 
    
   
    legend.position = "right",             
    legend.box = "vertical",               
    legend.margin = margin(l = 15, unit = "pt"), 
    legend.spacing.y = unit(20, "pt"),           
    
   
    legend.title = element_text(face = "bold", size = 11, color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(face = "bold", color = "black", size = 13), 
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b=15)),
    
    
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray95")
  )

pdf("WT-SMC.pdf", width = 8.5, height = 8)
print(p)
dev.off()

