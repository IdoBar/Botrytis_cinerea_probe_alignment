devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "util.R")
bioc_packs <- c("GenomicRanges", "rtracklayer", "GenomicFeatures", "Rsamtools", 
                "Gviz")
BiocManager::install(bioc_packs, update = FALSE)
`BiocManager::install("Gviz")

pacman::p_load(char = c("tidyverse","seqinr", bioc_packs ))


#### Read data ####
primer_table <- readxl::read_excel("./data/primer_info_marzia.xlsx") 
primer_table %>% filter(Target_species=="B. cinerea") %>% group_by(NCBI_Accession) %>% dplyr::slice(1) %>% 
  dplyr::select(Gene, NCBI_Accession) %>%   write_tsv("./data/primer_genes.txt")
# export as fasta for BLAST
util$write.fasta(primer_table[c("Primer", "Primer_sequence")], 
            filename = "data/primer_seqs.fasta")
# BLAST the file against the genome
# read blast results
blast_res <- read_tsv("data/primer_seqs.outfmt6", col_names = c("qseqid", "sseqid", "pident", "length",
"mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")) 


# 
#   dplyr::filter(!is.na(NCBI_Accession)) %>% 
#   dplyr::rename(Gene_coords=`Gene genomic coords (Me14)`) %>% 
#   separate(Gene_coords, 
#            c("Chr", "Gene_start", "Gene_end"), sep = "([\\:-])", remove = FALSE) %>% 
#   separate(`Primer coords`, c("Primer_start", "Primer_end")) %>% 
#   mutate(Amplicon= sub("[-/][RF]$", "", Primer)) %>% 
#   mutate_at(vars(ends_with("_end")), as.numeric) %>% 
#   mutate_at(vars(ends_with("_start")), as.numeric)

wgs_folder <- "../AGRF_CAGRF20074_HFLHKDRXX/AGRF_BT2_process_19_09_2019"
bam_files <- list.files(wgs_folder , pattern = ".bam$", full.names = TRUE)
# bam_files %>% walk(~indexBam(.x))
# indexBam(bam_files)

gff_file <- "data/genome_assemblies/GCA_000143535.4_ASM14353v4_genomic.gff"
options(ucscChromosomeNames=FALSE)
# import genome as fasta
seqs <- Biostrings::readDNAStringSet("data/genome_assemblies/GCA_000143535.4_ASM14353v4_genomic.fna")

# Import gff as GRanges
# gff <- import.gff(gff_file)

# Import gff as TxDb
txdb <- GenomicFeatures::makeTxDbFromGFF(gff_file, format="gff3", dataSource = "GCA_000143535.4",
                        organism = "Botrytis cinerea B05.10",
                        taxonomyId = 332648)

primer_blast <- primer_table %>% filter(grepl("cinerea", Target_species)) %>% 
  left_join(blast_res %>% dplyr::select(Chr=sseqid, Primer=qseqid))
accessions <- unique(primer_blast$NCBI_Accession)
for (acc in accessions){
  # acc <- accessions[1]
  reg <- primer_blast %>% dplyr::filter(NCBI_Accession==acc)
  loc <- unique(reg$Gene_coords)
  # st <- unique(reg$Gene_start) # + reg$Primer_start-5
  # en <- unique(reg$Gene_end) # + reg$Primer_end+5
  chr <- unique(reg$Chr)
  chr_name <- names(seqs)[grepl(unique(reg$Chr), names(seqs))]
  #### Verify primer sequence ####
  seq <- seqs[[chr_name]]
  
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
    seqs, genome = "Botrytis cinerea B05.10",
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
  
  reg_bam <- bam_files[grepl(glue::glue("cinerea.+{chr}"), bam_files)] 
  altrack_Bc <- AlignmentsTrack(
    bam_files[grepl(glue::glue("cinerea.+{chr}"), bam_files)], isPaired = TRUE,  
    showMismatches=TRUE, name = "B. cinerea reads",# col.mates = "deeppink",
    # min.height = 12, max.height = 20, coverageHeight = 0.15, size = 50, 
    reverseStacking=FALSE, min.height = 2.5, # coverageHeight = 0.15,
    type=c("coverage", "pileup"), start = st, end = en
  )
  
  altrack_Bf <- AlignmentsTrack(
    bam_files[grepl(glue::glue("fabae.+{chr}"), bam_files)], isPaired = TRUE,  
    showMismatches=TRUE, name = "B. fabae reads",# col.mates = "deeppink", 
    # min.height = 12, max.height = 20, coverageHeight = 0.15, size = 50, 
    reverseStacking=FALSE, min.height = 2.5, # coverageHeight = 0.15,
    type=c("coverage", "pileup"), start = st, end = en
  )
  # add highlighting 
  ht1 <- HighlightTrack(trackList=list(altrack_Bc, altrack_Bf), 
                        start = c(primers_pos$start[1], primers_pos$start[2]),
                        end = c(primers_pos$end[1],primers_pos$end[2]), 
                        chromosome=chr)
  LogMsg(glue::glue("Preparing plot for {reg$Gene[1]} on {chr}"))
  # plot 
  pdf(glue::glue("output/plots/Gviz_alignment_{reg$Gene[1]}_region_{chr}_{st}-{en}.pdf"),
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
