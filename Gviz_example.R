library("Gviz")
library("GenomicFeatures")
library(readr)

txdb <- makeTxDbFromGFF("at3g48187.gff3", format="gff3")
b <- GeneRegionTrack(txdb, col.line="darkgrey", name="Genome", transcriptAnnotation = "isoform", fill  ="salmon1")
plotTracks(b)


positionsdf <-  read_csv("at3g48187_rpg_len_7_35_+_positions_df.csv")

#TRYPSIN
trypsindf <- subset(positionsdf, enzyme=='Trypsin')
trypsin <- AnnotationTrack(start = trypsindf$gen_start_pos, end=trypsindf$gen_end_pos, 
                           chromosome = "chr3", 
                           strand = "+",
                           col.line = "white",
                           featureAnnotation=NULL,
                           group = trypsindf$isoform,
                           fill = "#F0E442",
                             name = "Trypsin")

#ARGC
argcdf <- subset(positionsdf, enzyme=='Arg-C')
argc <- AnnotationTrack(start = argcdf$gen_start_pos, end=argcdf$gen_end_pos, 
                        chromosome = "chr3", 
                        strand = "+",
                        col.line = "white",
                        featureAnnotation=NULL,
                        group = argcdf$isoform,
                        fill = "#009E73",
                        name = "Arg-C")

#ASPN
aspndf <- subset(positionsdf, enzyme=='Asp-N')
aspn <- AnnotationTrack(start = aspndf$gen_start_pos, end=aspndf$gen_end_pos, 
                        chromosome = "chr3", 
                        strand = "+",
                        col.line = "white",
                        featureAnnotation=NULL,
                        group = aspndf$isoform,
                        fill = "#56B4E9",
                        name = "Asp-N")

#CHYMOTRYPSIN
chymdf <- subset(positionsdf, enzyme=='Chymotrypsin-high')
chym <- AnnotationTrack(start = chymdf$gen_start_pos, end=chymdf$gen_end_pos, 
                        chromosome = "chr3", 
                        strand = "+",
                        col.line = "white",
                        featureAnnotation=NULL,
                        group = chymdf$isoform,
                        cex.title = 0.7,
                        fill = "#D55E00",
                        name = "Chym")

#GluC
glucdf <- subset(positionsdf, enzyme=='Glu-C')
gluc <- AnnotationTrack(start = glucdf$gen_start_pos, end=glucdf$gen_end_pos, 
                        chromosome = "chr3", 
                        strand = "+",
                        col.line = "white",
                        featureAnnotation=NULL,
                        group = glucdf$isoform,
                        fill = "#E69F00",
                        name = "Glu-C")

#LYSC
lyscdf <- subset(positionsdf, enzyme=='Lys-C')
lysc <- AnnotationTrack(start = lyscdf$gen_start_pos, end=lyscdf$gen_end_pos, 
                        chromosome = "chr3", 
                        strand = "+",
                        col.line = "white",
                        featureAnnotation=NULL,
                        group = lyscdf$isoform,
                        fill = "lightpink",
                        name = "Lys-C")

#LYSN
lysndf <- subset(positionsdf, enzyme=='Lys-N')
lysn <- AnnotationTrack(start = lysndf$gen_start_pos, end=lysndf$gen_end_pos, 
                        chromosome = "chr3", 
                        strand = "+",
                        col.line = "white",
                        featureAnnotation=NULL,
                        group = lysndf$isoform,
                        fill = "#0072B2",
                        name = "Lys-N")

plotTracks( c(b, argc, aspn, chym, gluc, lysc, lysn, trypsin), just.group="right", shape = "box", size=1, cex.group=1, cex.title=1, groupAnnotation = "group")
