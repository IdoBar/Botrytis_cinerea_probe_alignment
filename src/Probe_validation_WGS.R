devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "util.R")
BiocManager::install(c("GenomicRanges", "rtracklayer", "GenomicFeatures", "Rsamtools", 
                       "Gviz"), update = FALSE)
BiocManager::install("Gviz")
library(Gviz)
library(tidyverse)
pacman::p_load(tidyverse)


#### Read data ####
primer_table <- readxl::read_excel("./A.rabiei_primers/Primer_table.xlsx") %>% 
  dplyr::filter(!is.na(NCBI_Accession)) %>% 
  dplyr::rename(Gene_coords=`Gene genomic coords (Me14)`) %>% 
  separate(Gene_coords, 
           c("Chr", "Gene_start", "Gene_end"), sep = "([\\:-])", remove = FALSE) %>% 
  separate(`Primer coords`, c("Primer_start", "Primer_end")) %>% 
  mutate(Amplicon= sub("-[RF]$", "", Primer)) %>% 
  mutate_at(vars(ends_with("_end")), as.numeric) %>% 
  mutate_at(vars(ends_with("_start")), as.numeric)

bam_files <- list.files("./alignment_files/" , pattern = ".bam$", full.names = TRUE)

gff_file <- "Reference_genomes/Arab_me14.gff3"
options(ucscChromosomeNames=FALSE)
# import genpome as fasta
seqs <- Biostrings::readDNAStringSet("Reference_genomes/A.rabiei_me14.fasta")

# Import gff as GRanges
# gff <- import.gff(gff_file)

# Import gff as TxDb
txdb <- GenomicFeatures::makeTxDbFromGFF(gff_file, format="gff3", dataSource = "CCDM Curtin",
                        organism = "Ascochyta rabiei",
                        taxonomyId = 5454)


accessions <- unique(primer_table$NCBI_Accession)
for (acc in accessions){
  # acc <- accessions[1]
  reg <- primer_table %>% dplyr::filter(NCBI_Accession==acc)
  loc <- unique(reg$Gene_coords)
  # st <- unique(reg$Gene_start) # + reg$Primer_start-5
  # en <- unique(reg$Gene_end) # + reg$Primer_end+5
  chr <- unique(reg$Chr)
  #### Verify primer sequence ####
  seq <- seqs[[chr]]
  
  primers_pos <- str_locate(toString(seq),  reg$Primer_sequence) %>% as_tibble() 
  missing_pos <- which(rowSums(is.na(primers_pos))>0)
  primers_pos[missing_pos,] <- str_locate(toString(seq),  
              stringi::stri_reverse(chartr("ATGC","TACG", reg$Primer_sequence[missing_pos])))
  st <- min(primers_pos) - 20   
  en <- max(primers_pos) + 20
  #### Plot Regions ####
  # Test
  # ggplot() + geom_alignment(gff, type = "type")
  
  # Plot with Gviz
  # Add a sequence track
  strack <- SequenceTrack(
    seqs, genome = "A. rabiei Me14",
    chromosome = chr, from = st, to = en,
    cex=0.5, min.height =7, legend=TRUE
  )
  # add a genome axis
  gtrack <- GenomeAxisTrack()
  # add annotation
  grtrack <- GeneRegionTrack(
    txdb,
    chromosome = chr, start = st, end = en,
    showId = TRUE, stacking = "dense",
    name = "Gene Annotation"
  )
  # add alignment
  
  reg_bam <- bam_files[grepl(glue::glue('{sub(":", "_", loc, fixed = TRUE)}'), bam_files)]
  altrack <- AlignmentsTrack(
    reg_bam, isPaired = TRUE,  showMismatches=TRUE, # col.mates = "deeppink",
    # min.height = 12, max.height = 20, coverageHeight = 0.15, size = 50, 
    reverseStacking=FALSE, min.height = 2.5, # coverageHeight = 0.15,
    type=c("coverage", "pileup"), start = st, end = en
  )
  
  # add highlighting 
  ht1 <- HighlightTrack(trackList=list(altrack), 
                        start = c(primers_pos$start[1], primers_pos$start[2]),
                        end = c(primers_pos$end[1],primers_pos$end[2]), 
                        chromosome=chr)
  LogMsg(glue::glue("Preparing plot for {reg$Amplicon[1]} on {chr}"))
  # plot 
  pdf(glue::glue("plots/Gviz_alignment_{reg$Amplicon[1]}_region_{chr}_{st}-{en}.pdf"),
      width = 8, height = 6)
  plotTracks(
    list(gtrack, strack, ht1, grtrack), # 
    from = st, to = en
    # background.title="lightskyblue"
  )
  dev.off()
}

#### focus on single primer region ####
for (p in primer_table$Primer){
  # p <- primer_table$Primer[2]
  reg <- primer_table %>% dplyr::filter(Primer==p)
  loc <- unique(reg$Gene_coords)
  # st <- unique(reg$Gene_start) # + reg$Primer_start-5
  # en <- unique(reg$Gene_end) # + reg$Primer_end+5
  chr <- unique(reg$Chr)
  #### Verify primer sequence ####
  seq <- seqs[[chr]]
  
  primers_pos <- str_locate(toString(seq),  reg$Primer_sequence) %>% as_tibble() 
  missing_pos <- which(rowSums(is.na(primers_pos))>0)
  if (length(missing_pos)>0) {
    primers_pos[missing_pos,] <- str_locate(toString(seq),  
                                 stringi::stri_reverse(chartr("ATGC","TACG", 
                                                       reg$Primer_sequence[missing_pos])))
  }
  
  st <- min(primers_pos) - 20   
  en <- max(primers_pos) + 20
  #### Plot Regions ####
  # Test
  # ggplot() + geom_alignment(gff, type = "type")
  
  # Plot with Gviz
  # Add a sequence track
  strack <- SequenceTrack(
    seqs, genome = "A. rabiei Me14",
    chromosome = chr, from = st, to = en,
    cex=0.8, min.height =7, legend=TRUE
  )
  # add a genome axis
  gtrack <- GenomeAxisTrack()
  # add annotation
  grtrack <- GeneRegionTrack(
    txdb,
    chromosome = chr, start = st, end = en,
    showId = TRUE, stacking = "dense",
    name = "Gene Annotation"
  )
  # add alignment
  
  reg_bam <- bam_files[grepl(glue::glue('{sub(":", "_", loc, fixed = TRUE)}'), bam_files)]
  altrack <- AlignmentsTrack(
    reg_bam, isPaired = TRUE,  showMismatches=TRUE, # col.mates = "deeppink",
    # min.height = 12, max.height = 20, coverageHeight = 0.15, size = 50, 
    reverseStacking=FALSE, min.height = 2.5, # coverageHeight = 0.15,
    type=c("coverage", "pileup"), start = st, end = en
  )
  
  # add highlighting 
  ht1 <- HighlightTrack(trackList=list(altrack), 
                        start = c(primers_pos$start[1]),
                        end = c(primers_pos$end[1]), 
                        chromosome=chr)
  LogMsg(glue::glue("Preparing plot for {reg$Primer} on {chr}"))
  # plot 
  pdf(glue::glue("plots/Gviz_alignment_{reg$Primer}_region_{chr}_{st}-{en}.pdf"),
      width = 8, height = 6)
  plotTracks(
    list(gtrack, strack, ht1, grtrack), # 
    from = st, to = en
    # background.title="lightskyblue"
  )
  dev.off()
}

# check region with high polymorphism to compare

st <- 1871550# 1884621 # + reg$Primer_start-5
en <- 1871650# 1884675 # + reg$Primer_end+5
chr <- "Arab_Me14_ctg07"
#### Highly Polymorphic Region ####
# Test
# ggplot() + geom_alignment(gff, type = "type")

# Plot with Gviz
# Add a sequence track
strack <- SequenceTrack(
  seqs, genome = "A. rabiei Me14",
  chromosome = chr, from = st, to = en,
  cex=0.8, min.height =7
)
# add a genome axis
gtrack <- GenomeAxisTrack()
# add annotation
grtrack <- GeneRegionTrack(
  txdb,
  chromosome = chr, start = st, end = en,
  showId = TRUE, stacking = "dense",
  name = "Gene Annotation"
)
# add alignment
reg_bam <- bam_files[grepl(glue::glue("{chr}_1871289"), bam_files)]
altrack <- AlignmentsTrack(
  reg_bam, isPaired = TRUE,  showMismatches=TRUE, # col.mates = "deeppink",
  # min.height = 12, max.height = 20, coverageHeight = 0.15, size = 50, 
  reverseStacking=FALSE, min.height = 2.5, # coverageHeight = 0.15,
  type=c("coverage", "pileup"), stacking="squish", start = st, end = en, chromosome = chr
)

# add highlighting 
LogMsg(glue::glue("Preparing plot for polymorphic region on {chr}"))
# plot 
pdf(glue::glue("plots/Gviz_alignment_polymorphic_region_{chr}_{st}-{en}.pdf"), width = 8,
    height = 6)
plotTracks(
  list(gtrack, strack, altrack, grtrack), # 
  from = st, to = en
  # background.title="lightskyblue"
)
dev.off()

#### Mating-Type region ####
MT_region <- c(1667555, 1668666)
plot_region <- c(1666650, 1667750) # 850
reg <- primer_table %>% dplyr::filter(Gene=="Mating-Type")
loc <- unique(reg$Gene_coords)
# st <- unique(reg$Gene_start) # + reg$Primer_start-5
# en <- unique(reg$Gene_end) # + reg$Primer_end+5
chr <- unique(reg$Chr)
#### Verify primer sequence location ####
seq <- seqs[[chr]]

primers_pos <- str_locate(toString(seq),  reg$Primer_sequence) %>% as_tibble() 
missing_pos <- which(rowSums(is.na(primers_pos))>0)
primers_pos[missing_pos,] <- str_locate(toString(seq),  
                                        stringi::stri_reverse(chartr("ATGC","TACG", 
                                               reg$Primer_sequence[missing_pos])))
reg[c("Primer_start", "Primer_end")]
st <- plot_region[1]
en <- plot_region[2]
# Test
# ggplot() + geom_alignment(gff, type = "type")

# Plot with Gviz
# set colours
pal_cols <- paletteer::paletteer_d("dutchmasters", "view_of_Delft")
bar_col <- pal_cols[1]
left_bar_col <- pal_cols[3]
left_bar_text <- paletteer::paletteer_d("dutchmasters", "pearl_earring")[5] # pal_cols[4]
# calculate GC and add as a data track
calc_GC <- function(DNAseq, i, window=100){
  # DNAseq=seq
  # i=10
  start_win <- max(1, i-window/2)
  end_win <- min(length(DNAseq), i+window/2)
  GC <- Biostrings::letterFrequency(substring(DNAseq, start_win, end_win), letters = "GC", 
                        as.prob = TRUE)
}

plot_seq <- seq[st:en]
win_width=20
GC_win <- Biostrings::letterFrequencyInSlidingView(plot_seq, win_width, letters = "GC", 
                                                   as.prob = TRUE)

# GC_win_tbl <- tibble(start=seq_along(GC_win), end = seq_along(GC_win)+win_width-1)
# 
# GC_data <- tibble(GC=seq_along(plot_seq) %>% map_dbl(~calc_GC(plot_seq, i=.)))

GC_region <- IRanges(start=st-1+seq_along(GC_win),end=st-1+seq_along(GC_win)+win_width-1)
mcols(GC_region) <- data.frame(GC=GC_win)
GC_region <- GRanges(seqnames=chr,ranges = GC_region)
# Add a data track
dtrack <- DataTrack(GC_region, genome = "A. rabiei Me14", name = "GC",type="a",
                    background.title=left_bar_col, ylim = c(0,1),
                    fontcolor.title = left_bar_text) # , type="l"
# Add a sequence track
strack <- SequenceTrack(
  seqs, genome = "A. rabiei Me14",
  chromosome = chr, from = st, to = en,
  cex=0.8, min.height =7
)
# add a genome axis
gtrack <- GenomeAxisTrack(range=IRanges(start=st, end=en,names=chr),
                          labelPos="below", cex.id=1) # showId=TRUE, fill.range=bar_col,
# add annotation
grtrack <- GeneRegionTrack(
  txdb, transcriptAnnotation="symbol",
  chromosome = chr, start = st, end = en,
  showId = TRUE, stacking = "dense", background.title="transparent",
  name = "Gene Annotation"#, min.height =10
)
# add alignment
reg_bam <- bam_files[grepl(glue::glue('{sub(":", "_", loc, fixed = TRUE)}'), bam_files)]
altrack <- AlignmentsTrack(
  reg_bam, isPaired = TRUE,  showMismatches=TRUE, # col.mates = "deeppink",
  chromosome = chr, start = st, end = en,background.title=left_bar_col, 
  fontcolor.title = left_bar_text,
  col.axis=left_bar_text,
   
  # min.height = 12, max.height = 20, coverageHeight = 0.15, size = 50, 
  reverseStacking=FALSE, min.height = 2.5, # coverageHeight = 0.15,
  type=c("coverage", "pileup")
)

# add highlighting 
ht1 <- HighlightTrack(trackList=list(altrack), 
            start = c(primers_pos$start[1], primers_pos$start[2]),
            end = c(primers_pos$end[1],  primers_pos$end[2]), 
            # col=c("#E41A1C", "#377EB8"), 
            # fill=c("#FBB4AE", "#B3CDE3"),
                      # inBackground=FALSE,
                      chromosome=chr)
# ht2 <- HighlightTrack(trackList=list(altrack), col="royalblue",
#                       start = MT_region[1],
#                       end = MT_region[2], 
#                       chromosome=chr)
# displayPars(ht1) <- list(alpha.title=1, alpha=0.5)
# displayPars(ht2) <- list(alpha.title=1, alpha=0.5)
# ot <- OverlayTrack(trackList=list(ht1, ht2))
LogMsg(glue::glue("Preparing plot for {reg$Amplicon[1]} on {chr}"))
# plot 

pdf(glue::glue("plots/Gviz_alignment_GC_{reg$Amplicon[1]}_region_{chr}_{st}-{en}.pdf"),
    width = 8, height = 6)
plotTracks(
  list(gtrack, dtrack, strack, ht1, grtrack), # 
  from = plot_region[1], to = plot_region[2] 
)
dev.off()
