workprefix = 'otu_sharing'
options(prompt=paste0(workprefix, "> "), stringsAsFactors=FALSE)
reload = function(doit=FALSE) if (doit) source(paste0(workprefix, ".R"))


library(RColorBrewer)
library(Cairo)

dataprefix = "edgar_20160602_first_200_1.5"
pfile = paste0("otutab.", dataprefix, ".txt")
#ptree.all = "ttt.tre"
#ptree.OTUs.default = "sanger-plus-OTUs.muscle.out.default.no-sanger.degap.tre"
#ptree.OTUs.manual = "sanger-plus-OTUs.muscle.out.manual.20141214.no-sanger.degap.tre"

orig.otu.names = c("OTU", "Ice_TV", "Ice_CI", "Ol1_TV", "Ol1_CI", "Ol2_TV", "Ol2_CI", "Control")
new.otu.names = c("Ice_TV", "Ice_CI", "Ol1_TV", "Ol1_CI", "Ol2_TV", "Ol2_CI")
new.otu.longsites = c("Iceland", "Iceland", "Öland 1", "Öland 1", "Öland 2", "Öland 2")
new.otu.longspecies = c("Thamnolia", "Cetraria", "Thamnolia", "Cetraria", "Thamnolia", "Cetraria")
names(new.otu.longsites) = names(new.otu.longspecies) = new.otu.names


read.otu.matrix = function(file) {
    dat = read.delim(file, header=TRUE)
    names(dat) = orig.otu.names  # overwrites existing row names
    dat$OTU = gsub(paste0("_", dataprefix, "_"), " ", dat$OTU)
    rownames(dat) = dat$OTU
    dat$OTU = NULL  # drop OTU name column (duplicates row names)
    dat$Control = NULL
    ####
    dat
}

otu.sparklines = function(dat, do.pdf=FALSE) {
    at.midpoint = TRUE  # center boxes vertically at their midpoints
    dat = dat[nrow(dat):1, ]  # reverse rows
    # scale each sample by its sum, so transform to [0, 1]
    dat = scale(dat, center=FALSE, scale=(apply(dat, 2, sum)))
    this.xlim = c(0, ncol(dat))
    this.ylim = c(0, nrow(dat))
    xbase = 1:ncol(dat)
    ybase = 1:nrow(dat)
    yscale = 1.5 #1 #0.8

    if (do.pdf) CairoPDF(file=paste0("otu_sparklines_", dataprefix, ".pdf"), width=4.5, height=6)

    def.par = par(no.readonly=TRUE)
    layout(matrix(c(0, 2, 3, 1), 2, 2, byrow=TRUE),
           widths = c(0.6, 3), heights = c(0.5, 3))
    opa.1 = par(mar=c(0.0, 1.5, 0, 0.0), mgp=c(1.0, 0, 0), las=1)
    colfunc = colorRampPalette(c("brown2","darkblue"))
    cols = colfunc(nrow(dat))
    plot.new()
    plot.window(xlim=this.xlim, ylim=this.ylim)

    clustering.method = "canberra"

    sites.clust = as.dendrogram(hclust(dist(t(dat), method=clustering.method)))
    sites.clust = rev(reorder(sites.clust, c(6:1), mean))
    # reorder columns
    reordered.sites = c("Ice_TV", "Ice_CI", "Ol1_CI", "Ol2_CI", "Ol1_TV", "Ol2_TV")
    dat = dat[, reordered.sites]
    new.otu.longsites = new.otu.longsites[reordered.sites]
    new.otu.longspecies = new.otu.longspecies[reordered.sites]

    otus.clust = as.dendrogram(hclust(dist(dat, method=clustering.method)))
    otus.clust = rev(otus.clust)
    # reorder rows
    reordered.rows = labels(otus.clust)
    dat = dat[reordered.rows, ]

    attributes(dat)$cluster.sites = sites.clust
    attributes(dat)$cluster.otus = otus.clust

    for (r in ybase) {
        d = data.frame(x=xbase, y=as.vector(t(dat[r,])))
        d = subset(d, !is.na(y) & y > 0)
        if (nrow(d) == 0)
            next
        x0 = d$x - 1; x1 = d$x
        y0 = r - 1;   y1 = y0 + (d$y^0.333)
        if (at.midpoint) {
            ymid = y0 + 0.5
            dy2 = (y1 - y0)/2
            y0 = ymid - dy2
            y1 = ymid + dy2
        }
        rect(x0, y0, x1, y1, col=cols[r], border=NA)
    }

    y0 = max(ybase) + 3.0
    text(x=xbase - 0.5, y=y0, labels=new.otu.longsites, pos=1, cex=0.8, xpd=NA)
    text(x=xbase - 0.5, y=y0 - 1, labels=new.otu.longspecies, font=3, pos=1, cex=0.8, xpd=NA)

    # OTU names at left side, and guide lines below names
    x0 = -0.8
    text(x=x0, y=ybase - 0.5, labels=rownames(dat), pos=4, cex=0.6, xpd=NA)
    y = ybase - 1

    # clustering; euclidean too dependent on frequency and clusters rare OTUs
    # better to use 'canberra' (better) or 'binary' which are more presence/absense
    # 
    # how do we represent phylogeny in this?  or do we?
    # maybe cluster by site in this figure, then leave OTUs as they are in the table
    # then a separate figure which plots OTU phylogeny together with OTU clustering based on sharing
    #     the dendextend packages can help
    #

    # plot at top
    opa2 = par(mar=c(1.5, 1.5, 0.0, 0.0))
    plot(attr(dat, 'cluster.sites'), leaflab="none", yaxt="n", xpd=NA)

    # plot at left
    opa2 = par(mar=c(0.00, 0.0, 0.0, 0.0))
    plot(attr(dat, 'cluster.otus'), horiz=TRUE, leaflab="none", yaxt="n", xpd=NA, cex=0.4)

    if (do.pdf) dev.off()

    invisible(dat)
}

