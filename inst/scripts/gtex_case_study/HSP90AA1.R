rm(list=ls())
path="C:/Document/Serieux/Travail/Data_analysis_and_papers/GTEX_analysis_Fsusie"
source(paste0(path,"/code/plot_all_effect_log.R"))
source(paste0(path,"/code/plot_log.R"))

load(paste0(path,"/output/local/ENSG00000080824.csv.gz.RData")) 

library(ggplot2)
library(fsusieR)

plot_susiF_pip(out$res, pos_SNP = as.numeric(out$info_SNP$POS))+
  
  geom_rect(aes(xmin = out$locus[1], 
                xmax =out$locus[2],
                ymin = -0.01,
                ymax = 1.02), 
            alpha = 0.0, color = "red")+
  xlim(c(out$locus[1] -150000,out$locus[2] +150000 )  ) 





folder_path=  paste0(getwd(),
                     "/plot/HSP90AA1/"
)
file_path <- file.path(folder_path, "HSP90AA1_pip.pdf")
pdf(file_path, width =  11.69, height =8.27)  # A4 in inches

plot_susiF_pip(out$res, pos_SNP = as.numeric(out$info_SNP$POS))+
  
  geom_rect(aes(xmin = out$locus[1], 
                xmax =out$locus[2],
                ymin = -0.01,
                ymax = 1.02), 
            alpha = 0.0, color = "red")+
  xlim(c(out$locus[1] -150000,out$locus[2] +150000 )  ) +
  theme(
    axis.text = element_text(size = 14),       # Increase axis text size
    axis.title = element_text(size = 16, face = "bold"), # Increase axis labels size
    legend.text = element_text(size = 14),     # Increase legend text size
    legend.title = element_text(size = 16, face = "bold"), # Increase legend title size
    egend.key.size = unit(1.5, "cm"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5) # Increase title size
  )

dev.off()

obj= out$res
Y=log1p(as.matrix(out$Y/out$size_factor_local) )
X=as.matrix(out$X) 
 


m1 <-susieR::susie(X=X,
                   y=rowSums(Y),
                   L=5
)
m1$sets
 

PCA <- svd(Y)
m2 <-susieR::susie(X=X,
                   y=PCA$u[,1],
                   L=5
)

m2$sets




out1=out
obj1=out$res
obj=out$res
chr= paste0("chr",out1$chr)
pos0 =out1$locus[1]
pos1=out1$locus[2]

snp_info=out1$info_SNP
cs = 1
log1p_count=TRUE
data_splice=NULL
plot_cred_band=TRUE
type_data="p"

# Extract the relevant genes and exons in the specified region
region_genes <- genes(txdb,columns = c("tx_id","gene_id"))

# Subset the genes and exons to the region of interest.
region_genes <- subsetByOverlaps(region_genes,
                                 GRanges(seqnames = chr,
                                         ranges = IRanges(pos0,pos1)))

# Generate a sequence of positions with a length of 1,024.
positions <- seq(pos0,pos1,length.out = 1024)

markers <- obj$cs[[cs]]
j       <- which.max(obj$pip[markers])
marker  <- markers[j]
x       <- X[,marker]





read_counts <- rbind(colMeans(log1p(out$Y[x == 0,])),
                     colMeans(log1p(out$Y[x == 1,])) )

### CS1 -----

uni_res= univariate_functional_regression(Y= fsusieR:::colScale(Y  )   ,
                                          X=as.matrix(X[ ,obj$cs[[cs]][1]],
                                                      ncol=1) ,
                                          method="TI",alpha = 0.1
)

rest =smash_regression(obj, Y=   fsusieR:::colScale(Y  , scale = FALSE)  ,
                     X=  X)#fsusieR:::colScale(X    ))

plot(rest$fitted_func[[1]])

obs_diff=read_counts[2,]-read_counts[1,]
plot(smashr::smash(obs_diff))

effect=rbind (uni_res$effect_estimate ,
              uni_res$cred_band ,
              rep(0, length(obj$fitted_func[[cs]])),
              obs_diff)



effect= rbind(rest$fitted_func[[1]]/2,
              rest$cred_band[[1]] /2,
              rep(0, length(obj$fitted_func[[cs]])) )
t_effect=  smashr::smash(read_counts[2,]-read_counts[1,], post.var=TRUE,
                         sigma=sqrt(var(read_counts[2,]-read_counts[1,])))

cred_band = rbind( t_effect$mu.est+1.98*sqrt(t_effect$mu.est.var),
                   t_effect$mu.est-1.98*sqrt(t_effect$mu.est.var))
effect=rbind (t_effect$mu.est,
              cred_band ,
              rep(0, length(obj$fitted_func[[cs]])) )

obs_effect=rbind (read_counts[2,]-read_counts[1,],
                  rep(0, length(obj$fitted_func[[cs]])) )

ylim= c( min( c(effect,obs_effect)), max(c(  effect,obs_effect)))
group_cred= c(1:3,0 )
group_colors <- c( "black"   ,"royalblue","royalblue","royalblue" )


# Create a "data track" to show the CS effect.
cex <- 1

group_lwd= c(1,2,1,1,1)
group_lty= c(1,1,2,2,1)
effect_track <-
  DataTrack(range = GRanges(seqnames = chr,
                            ranges = IRanges(start = positions,
                                             end = positions + 1)),
            data = effect, genome = "hg38",
            ylim=ylim,
            groups= group_cred,
            lwd=group_lwd,
            lty=group_lty,
            name = paste("Effect CS",cs),type = "l",col = group_colors,
            track.margin = 0.05,cex.title = cex,cex.axis = cex,
            col.axis = "black",col.title = "black",
            fontface = "plain",background.title = "white",
            fontface.title = 1,,
            legend = FALSE )


plotTracks(effect_track)
# Create another "data track" to show the read counts.


obs_effect=rbind (read_counts[2,]-read_counts[1,],
              rep(0, length(obj$fitted_func[[cs]])) )
 


obs_effect_track <-
  DataTrack(range = GRanges(seqnames = chr,
                            ranges = IRanges(start = positions,
                                             end = positions + 1)),
            data =read_counts[2,]-read_counts[1,], genome = "hg38",
            
            ylim=ylim,
            
            name ="Observed difference ",type = c("p"  ),col = c("royalblue"),
            track.margin = 0.05,cex.title = cex,cex.axis = cex,
            col.axis = "black",col.title = "black",
            fontface = "plain",background.title = "white",
            fontface.title = 1,,
            legend = FALSE )
tt=read_counts[2,]-read_counts[1,]

obs_effect_track2 <-
  DataTrack(range = GRanges(seqnames = chr,
                            ranges = IRanges(start = positions,
                                             end = positions + 1)),
            data =  rep(0, length(obj$fitted_func[[cs]])), genome = "hg38",
            
            ylim= ylim,
            name ="Observed difference ",type = c("l"  ),col = c("black" ),
            track.margin = 0.05,cex.title = cex,cex.axis = cex,
            col.axis = "black",col.title = "black",
            fontface = "plain",background.title = "white",
            fontface.title = 1,,
            legend = FALSE )

obs_effect_track =OverlayTrack(trackList = list(obs_effect_track,
                                                obs_effect_track2),background.title = "white")

plotTracks(obs_effect_track)

n0  <- sum(x == 0)
n1  <- sum(x == 1)
n2  <- sum(x == 2)
id  <- "rs376197939" #snp_info[marker,"ID"]
ref <- "AA"#snp_info[marker,"REF"]
alt <- "AT" #snp_info[marker,"ALT"]


groups <- c(sprintf("\t    %s %s (n = %d)",id, ref,n0),
            sprintf("\t    %s %s (n = %d)",id, alt,n1) )
geno_colors <- c("navyblue","turquoise" )



groups <- factor(groups,rev(groups))
geno_colors <- rev(geno_colors)


lab_y =ifelse(log1p_count, "Avg. log1p count","Avg. count")


data_track <-DataTrack(range = GRanges(seqnames = chr,
                                       ranges = IRanges(start = positions,
                                                        end = positions + 1)),
                       data = read_counts,genome = "hg38",
                       groups = groups,
                       name = lab_y , type = type_data,  col = geno_colors ,track.margin = 0.05,cex.title = cex,cex.axis = cex,
                       
                       col.axis = "black",col.title = "black",
                       fontface = "plain",background.title = "white", 
                       fontface.title = 1, cex= .6,cex.legend = 1.1)


plotTracks(data_track)

# Create an "ideogram" track.
ideo_track <- IdeogramTrack(genome = "hg38",chromosome = chr)

# Create a "genome axis" track.
genome_track <- GenomeAxisTrack(col.axis = "black",col.title = "black")

# Create a "gene region" track.
gene_track <- GeneRegionTrack(txdb,genome = "hg38",chromosome = chr,
                              pos0 = pos0,pos1 = pos1,name = "",
                              showId = TRUE,geneSymbol = TRUE,
                              col.axis = "black",col.title = "black",
                              transcriptAnnotation = "symbol",
                              rotation.title = 0,cex.title = 2,
                              col = "salmon",fill = "salmon",
                              background.title = "white")

# Map gene IDs to gene symbols.
gene_ids <- unique(unlist(region_genes$gene_id))

# Map to gene symbols using org.Hs.eg.db
gene_symbols <- AnnotationDbi::select(org.Hs.eg.db,keys = gene_ids,
                                      columns = "SYMBOL",
                                      keytype = "ENTREZID")
n <- nrow(gene_symbols)
if (n > 0) {
  for (i in 1:n) {
    j <- which(gene_track@range@elementMetadata@listData$gene ==
                 gene_symbols$ENTREZID[i])
    gene_track@range@elementMetadata@listData$id[j]     <- gene_symbols$SYMBOL[i]
    gene_track@range@elementMetadata@listData$symbol[j] <- gene_symbols$SYMBOL[i]
  }
}


# Combine all tracks into a single plot.
tracks <- c(ideo_track,
            genome_track,
            effect_track,
            obs_effect_track ,
            data_track,
            gene_track)


plotTracks(tracks,from = pos0,to = pos1,sizes = c(1,1.75,1.75,2,4,2)) 

folder_path=  paste0(getwd(),
                     "/plot/HSP90AA1/"
)
file_path <- file.path(folder_path, "HSP90AA1_cs1.pdf")
pdf(file_path, width =16.69, height = 8.27 )  # A4 in inches



plotTracks(tracks,from = pos0-200,
           to = pos1+200,
           sizes = c(1,1 ,4,4,4,2),
           cex.main=1.2, cex.title = 1.
) 
dev.off()







## CS 2 -----
cs=2
Y=log1p(as.matrix(out$Y/out$size_factor_local) )
X=as.matrix(out$X)

markers <- obj$cs[[cs]]
j       <- which.max(obj$pip[markers])
marker  <- markers[j]
x       <- X[,marker]
n0  <- sum(x == 0)
n1  <- sum(x == 1)
n2  <- sum(x == 2)
id  <- "rs1400511103" #snp_info[marker,"ID"]
ref <- "CC"#snp_info[marker,"REF"]
alt <- "CT" #snp_info[marker,"ALT"]


groups <- c(sprintf("    %s %s (n = %d)",id, ref,n0),
            sprintf("    %s %s (n = %d)",id, alt,n1) )
geno_colors <- c("lightgreen","darkgreen" )


read_counts <- rbind(colMeans(log1p(out$Y[x == 0,])),
                     colMeans(log1p(out$Y[x == 1,])) )

which(x==2
)# removing only individual with 2 variant



#effect=rbind (tt$fitted_func[[2]] ,
#              tt$cred_band[[2]] ,
#              rep(0, length(obj$fitted_func[[cs]])))

uni_res= univariate_functional_regression(Y= Y  ,
                                          X=as.matrix(X[,obj$cs[[2]][1]],
                                                      ncol=1),
                                          method="TI",alpha = 0.1
)


effect=rbind (uni_res$effect_estimate ,
              uni_res$cred_band ,
              rep(0, length(obj$fitted_func[[cs]])))
effect=rbind (rest$fitted_func[[2]]/2,
              rest$cred_band[[2]]/2 ,
              rep(0, length(obj$fitted_func[[cs]])))
t_effect=  smashr::smash(read_counts[2,]-read_counts[1,], post.var=TRUE,
                         sigma=sqrt(var(read_counts[2,]-read_counts[1,])))

cred_band = rbind( t_effect$mu.est+1.98*sqrt(t_effect$mu.est.var),
                   t_effect$mu.est-1.98*sqrt(t_effect$mu.est.var))
effect=rbind (t_effect$mu.est,
              cred_band ,
              rep(0, length(obj$fitted_func[[cs]])) )
obs_effect=rbind (read_counts[2,]-read_counts[1,],
                  rep(0, length(obj$fitted_func[[cs]])) )

ylim= c( min( c(effect,obs_effect)), max(c(  effect,obs_effect)))
group_cred= c(1:3,0)
group_colors <- c("black" ,"green4","green4","green4" )




# Create a "data track" to show the CS effect.
cex <- 1

group_lwd= c(1,2,1,1)
group_lty= c(1,1,2,2)
effect_track <-
  DataTrack(range = GRanges(seqnames = chr,
                            ranges = IRanges(start = positions,
                                             end = positions + 1)),
            data = effect, genome = "hg38",
            ylim=ylim,
            groups= group_cred,
            lwd=group_lwd,
            lty=group_lty,
            name = paste("Effect CS",cs),type = "l",col = group_colors,
            track.margin = 0.05,cex.title = cex,cex.axis = cex,
            col.axis = "black",col.title = "black",
            fontface = "plain",background.title = "white",
            fontface.title = 1,,
            legend = FALSE )


plotTracks(effect_track)
# Create another "data track" to show the read counts.


groups <- factor(groups,rev(groups))
# geno_colors <- rev(geno_colors)



lab_y =ifelse(log1p_count, "Avg. log1p count","Avg. count")

data_track <-DataTrack(range = GRanges(seqnames = chr,
                                       ranges = IRanges(start = positions,
                                                        end = positions + 1)),
                       data = read_counts,genome = "hg38",
                       groups = groups,
                       name = lab_y , type = type_data,  col = geno_colors ,
                       track.margin = 0.05,cex.title = cex,cex.axis = cex,
                       
                       col.axis = "black",col.title = "black",
                       fontface = "plain",background.title = "white", 
                       fontface.title = 1, cex= .6,cex.legend = 1.1)


plotTracks(data_track)





obs_effect=rbind (read_counts[2,]-read_counts[1,],
                  rep(0, length(obj$fitted_func[[cs]])) )
obs_effect_track <-
  DataTrack(range = GRanges(seqnames = chr,
                            ranges = IRanges(start = positions,
                                             end = positions + 1)),
            data =read_counts[2,]-read_counts[1,], genome = "hg38",
            ylim=ylim,
          
           
            name ="Observed difference ",type = c("p"  ),col = c("green4" ),
            track.margin = 0.05,cex.title = cex,cex.axis = cex,
            col.axis = "black",col.title = "black",
            fontface = "plain",background.title = "white",
            fontface.title = 1,,
            legend = FALSE )
tt=read_counts[2,]-read_counts[1,]

obs_effect_track2 <-
  DataTrack(range = GRanges(seqnames = chr,
                            ranges = IRanges(start = positions,
                                             end = positions + 1)),
            data =  rep(0, length(obj$fitted_func[[cs]])), genome = "hg38",
            
             ylim=ylim,
            name ="Observed difference ",type = c("l"  ),col = c("black" ),
            track.margin = 0.05,cex.title = cex,cex.axis = cex,
            col.axis = "black",col.title = "black",
            fontface = "plain",background.title = "white",
            fontface.title = 1,,
            legend = FALSE )

obs_effect_track =OverlayTrack(trackList = list(obs_effect_track,
                                                obs_effect_track2),,background.title = "white")
plotTracks(obs_effect_track )
# Combine all tracks into a single plot.
tracks <- c(ideo_track,
            genome_track,
            effect_track,
            obs_effect_track ,
            data_track,
            gene_track)


plotTracks(tracks,from = pos0,to = pos1,sizes = c(1,1.75,1.75,2,4,2)) 

folder_path=  paste0(getwd(),
                     "/plot/HSP90AA1/"
)
file_path <- file.path(folder_path, "HSP90AA1_cs2.pdf")
pdf(file_path, width =16.69, height = 8.27 )  # A4 in inches



plotTracks(tracks,from = pos0-200,
           to = pos1+200,
           sizes = c(1,1 ,4,4,4,2),
           cex.main=1.2, cex.title = 1.
) 
dev.off()


