# GO enrichment analysis
# Moisès Bernabeu
# Barcelona, May 2023
# Partial code from Saioa Manzano-Morales

library(topGO)
library(ggplot2)
library(dplyr)
library(stringr)
library(see)

theme_set(theme_bw())

go_enrichment <- function(go_data, genes) {
  gn_list <- factor(as.integer(gn_nms %in% genes))
  names(gn_list) <- gn_nms
  
  ontologies <- c("MF", "BP", "CC")
  res <- c()
  for (ont in ontologies) {
    # GOdata object
    GOdata <- new("topGOdata", ontology = ont, 
                  allGenes = gn_list, annot = annFUN.gene2GO, 
                  gene2GO = go_data)
    
    # Fisher's exact test with weighted algorithm
    test <- runTest(GOdata, algorithm = 'weight', statistic = 'fisher')

    # Result table
    # Note that we're keeping all results so that the FDR correction is precise
    allRes <- GenTable(GOdata, raw.p.value = test, topNodes = length(test@score), numChar = 120)
    allRes <- data.frame(allRes, 'Ontology' = ont)
    allRes$raw.p.value <- as.numeric(allRes$raw.p.value)
    res <- rbind(res, allRes)
  }
  
  return(res)
}

plot_go_enrichment <- function(res, odir, organism, label) {
  i <- 1
  for (ont in unique(res$Ontology)) {
    sub_df <- res[which(res$Ontology == ont), ]
    
    terms <- length(sub_df$Term)
    
    maxlen <- c()
    for (j in sub_df$Term) {maxlen <- c(maxlen, str_length(j))}
    maxlen <- max(maxlen)
    
    plot_wid <- 9/107 * maxlen
    plot_hei <- 8/53 * terms
    
    p <- ggplot(sub_df,
                aes(y = reorder(Term, -log10(raw.p.value)),
                    x = -log10(raw.p.value))) +
      geom_point(colour = okabeito_colors(i)) +
      geom_segment(aes(x = 0, xend = -log10(raw.p.value), y = Term, yend = Term),
                   colour = okabeito_colors(i)) +
      ylab('GO term description') +
      xlab('-log10(p-value)') +
      labs(title = ont)
    pdf(paste(odir, '/go_enrichment', organism, label, ont, '.pdf', sep = ''),
        width = plot_wid, height = plot_hei)
    print(p)
    dev.off()
    i <- i + 1
  }

  p <- ggplot(res,
              aes(y = reorder(Term, -log10(raw.p.value)),
                  x = -log10(raw.p.value), colour = Ontology)) +
    geom_point() +
    geom_segment(aes(x = 0, xend = -log10(raw.p.value), y = Term, yend = Term)) +
    ylab('GO term description') +
    xlab('-log10(p-value)') +
    scale_color_okabeito() +
    facet_grid(Ontology~., scale = 'free', space = 'free')
  
  print(p)
}

load('../data/sp2sp_dat.RData')

go_data <- readMappings('../data/seed_go/seed_goterms_ids.tsv', sep="\t", IDsep = ';')
gn_nms <- names(go_data)
names(gn_nms) <- gn_nms

ggplot(datc[which(datc$sp_to == 'PAPAN'), ], aes(ndist)) +
  geom_density() +
  geom_vline(xintercept = 1.3, lty = 4, col = 'darkorange3') +
  xlim(0, 5)

papan_short <- datc[which(datc$sp_to == 'PAPAN' & datc$ndist <= 1.3), 1]

papan_short_res <- go_enrichment(go_data, papan_short)
papan_short_res_p <- papan_short_res[which(papan_short_res$raw.p.value <= 0.01), ]

plot_go_enrichment(papan_short_res_p, '../outputs/', 'PAPAN', 'short')

ggplot(papan_short_res_p[which(papan_short_res_p$Annotated > 20), ],
       aes(y = reorder(Term, Significant / Annotated),
           x = Significant / Annotated, colour = raw.p.value)) +
  geom_point() +
  geom_segment(aes(x = 0, xend = Significant / Annotated, y = Term, yend = Term)) +
  ylab('GO term description') +
  xlab('Significant / Annotated') +
  facet_grid(Ontology~., scale = 'free', space = 'free') +
  labs(title = 'PAPAN Short branches')

papan_long <- datc[which(datc$sp_to == 'PAPAN' & datc$ndist > 1.3), 1]

papan_long_res <- go_enrichment(go_data, papan_long)
papan_long_res_p <- papan_long_res[which(papan_long_res$raw.p.value <= 0.01), ]

plot_go_enrichment(papan_long_res_p, '../outputs/', 'PAPAN', 'long')

ggplot(papan_long_res_p[which(papan_long_res_p$Annotated > 20), ],
       aes(y = reorder(Term, Significant / Annotated),
           x = Significant / Annotated, colour = raw.p.value)) +
  geom_point() +
  geom_segment(aes(x = 0, xend = Significant / Annotated, y = Term, yend = Term)) +
  ylab('GO term description') +
  xlab('Significant / Annotated') +
  facet_grid(Ontology~., scale = 'free', space = 'free') +
  labs(title = 'PAPAN Long branches')

# MACMU ----
ggplot(datc[which(datc$sp_to == 'MACMU'), ], aes(ndist)) +
  geom_density() +
  geom_vline(xintercept = 1.27, lty = 4, col = 'darkorange3') +
  xlim(0, 5)

macmu_short <- datc[which(datc$sp_to == 'MACMU' & datc$ndist <= 1.27), 1]

macmu_short_res <- go_enrichment(go_data, macmu_short)
macmu_short_res_p <- macmu_short_res[which(macmu_short_res$raw.p.value <= 0.01), ]

plot_go_enrichment(macmu_short_res_p, '../outputs/', 'MACMU', 'short')

ggplot(macmu_short_res_p[which(macmu_short_res_p$Annotated > 20), ],
       aes(y = reorder(Term, Significant / Annotated),
           x = Significant / Annotated, colour = raw.p.value)) +
  geom_point() +
  geom_segment(aes(x = 0, xend = Significant / Annotated, y = Term, yend = Term)) +
  ylab('GO term description') +
  xlab('Significant / Annotated') +
  facet_grid(Ontology~., scale = 'free', space = 'free') +
  labs(title = 'MACMU short branches')

macmu_long <- datc[which(datc$sp_to == 'MACMU' & datc$ndist > 1.3), 1]

macmu_long_res <- go_enrichment(go_data, macmu_long)
macmu_long_res_p <- macmu_long_res[which(macmu_long_res$raw.p.value <= 0.01), ]

plot_go_enrichment(macmu_long_res_p, '../outputs/', 'MACMU', 'long')

ggplot(macmu_long_res_p[which(macmu_long_res_p$Annotated > 20), ],
       aes(y = reorder(Term, Significant / Annotated),
           x = Significant / Annotated, colour = raw.p.value)) +
  geom_point() +
  geom_segment(aes(x = 0, xend = Significant / Annotated, y = Term, yend = Term)) +
  ylab('GO term description') +
  xlab('Significant / Annotated') +
  facet_grid(Ontology~., scale = 'free', space = 'free') +
  labs(title = 'MACMU long branches')
