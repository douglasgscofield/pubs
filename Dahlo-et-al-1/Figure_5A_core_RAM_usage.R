# Plots Figure 5A from Dahlö et al., core and RAM usage histograms.  Uses data from "Figure_5A_data_core_RAM_usage.RData"

options(prompt='Figure_5A> ', scipen=10)
reload <- function(doit=FALSE) if (doit) source('Figure_5A_core_RAM_usage.R')

library('stringr')
library('RColorBrewer')


read.data <- function(f='Figure_5A_data.RData') {
    load(f)
    hd.ngs_cpu <<- hd.ngs_cpu
    hd.oth_cpu <<- hd.oth_cpu
    hd.ngs_mem <<- hd.ngs_mem
    hd.oth_mem <<- hd.oth_mem
    cpu_log <<- cpu_log
    mem_log <<- mem_log
    weighted_data <<- weighted_data
    'read into environment:  hd.ngs_cpu  hd.oth_cpu  hd.ngs_mem  hd.oth_mem  cpu_log  mem_log  weighted_data'
}


# axis creation

do.plot <- function(do.pdf=FALSE, y.lines=TRUE, scale.counts=10^5) {

    n.ngs = sum(hd.ngs_cpu$counts)
    n.oth = sum(hd.oth_cpu$counts)

    if (is.numeric(scale.counts)) {
        # the y-axis scale becomes 'per scale.counts jobs'
        do.scale.counts = function(cnts) ((cnts / sum(cnts)) * scale.counts)
        hd.ngs_cpu$counts = do.scale.counts(hd.ngs_cpu$counts)
        hd.oth_cpu$counts = do.scale.counts(hd.oth_cpu$counts)
        hd.ngs_mem$counts = do.scale.counts(hd.ngs_mem$counts)
        hd.oth_mem$counts = do.scale.counts(hd.oth_mem$counts)
    }

    # log the counts
    hd.ngs_cpu$counts <- log10(hd.ngs_cpu$counts)
    hd.oth_cpu$counts <- log10(hd.oth_cpu$counts)
    hd.ngs_mem$counts <- log10(hd.ngs_mem$counts)
    hd.oth_mem$counts <- log10(hd.oth_mem$counts)

    xy = c(4.33, 2.55)
    r = 1200
    if (do.pdf) pdf(paste0('usage_hists', ifelse(is.numeric(scale.counts), paste0('_scaled-',scale.counts), ''), '.pdf'), xy[1], xy[2], pointsize=9)

    # plot options
    par(mfrow=c(2,2), las=1, mar=c(1, 1, 1, 0.5), oma=c(2, 2, 0.1, 0), mgp=c(1.5, 0.3, 0), tcl=-0.2)


    xlim.cpu = c(floor(min(cpu_log)), ceiling(max(cpu_log)))
    xlim.mem = c(floor(min(mem_log)), ceiling(max(mem_log)))
    ylim.cpu = c(0, max(0, hd.ngs_cpu$counts, hd.oth_cpu$counts, log10(scale.counts), na.rm=TRUE))
    ylim.mem = c(0, max(0, hd.ngs_mem$counts, hd.oth_mem$counts, log10(scale.counts), na.rm=TRUE))


    yat = seq(ylim.mem[1], ylim.mem[2], by=1)
    ylab = do.call(expression, lapply(yat, function(x) bquote(10^{.(x)})))

    xat_cpu = seq(floor(min(cpu_log)), ceiling(max(cpu_log)))
    xat_mem = seq(floor(min(mem_log)), ceiling(max(mem_log)))

    xlab_cpu = as.integer(2^seq(floor(min(cpu_log)), ceiling(max(cpu_log))))
    xlab_mem = 2^seq(floor(min(mem_log)), ceiling(max(mem_log)))/1024


    #ngs.col = brewer.pal(3,name='BuGn')[3]
    #oth.col = brewer.pal(3,name='YlOrBr')[3]
    oth.col = brewer.pal(9, "Reds")[7]
    ngs.col = brewer.pal(9, "Greens")[7]
    hist.lwd=0.75
    hist.lend=0

    my.plot.histogram(hd.ngs_cpu, xlim=xlim.cpu, ylim=ylim.cpu, col=ngs.col, border=ngs.col, lwd=hist.lwd, lend=hist.lend, xlab="", ylab="", axes=FALSE, main="")
    if (y.lines) segments(xlim.cpu[1]-1, yat, xlim.cpu[2], yat, lwd=0.5, col='grey', lty='dotted')
    points(log(16, base=2), yat[length(yat)]-0.5, pch=25, cex=0.8, bg='black') 
    axis(2, at=yat, labels=ylab) # y axis
    axis(1, at=xat_cpu, labels=FALSE, mgp=c(1.5, 0.4, 0.1))
    mtext("NGS projects", side=3, line=-0.2, font=2)

    my.plot.histogram(hd.ngs_mem, xlim=xlim.mem, ylim=ylim.mem, col=ngs.col, border=ngs.col, lwd=hist.lwd, lend=hist.lend, xlab="", ylab="", axes=FALSE, main="")
    if (y.lines) segments(xlim.mem[1]-1, yat, xlim.mem[2], yat, lwd=0.5, col='grey', lty='dotted')
    axis(2, at=yat, labels=FALSE)
    axis(1, at=xat_mem, labels=FALSE, mgp=c(1.5, 0.4, 0.1))
    #mtext(bquote(italic(N)[total]==.(n.ngs)*','~italic(N)[scaled]==10^.(log10(scale.counts))), side=3, adj=1, cex=0.8, line=-0.2)
    mtext(bquote(italic(N)[total]==.(n.ngs)~'(unscaled)'), side=3, adj=1, cex=0.8, line=-0.2)



    my.plot.histogram(hd.oth_cpu, xlim=xlim.cpu, ylim=ylim.cpu, col=oth.col, border=oth.col, lwd=hist.lwd, lend=hist.lend, xlab="", ylab="", axes=FALSE, main="")
    if (y.lines) segments(xlim.cpu[1]-1, yat, xlim.cpu[2], yat, lwd=0.5, col='grey', lty='dotted')
    #points(log(16, base=2), yat[length(yat)]-0.5, pch=25, cex=0.8, bg='black') 
    axis(2, at=yat, labels=ylab) # y axis
    axis(1, at=xat_cpu, labels=rep("", length(xat_cpu)), mgp=c(1.5, 0.4, 0.1)) # x axis
    text(xat_cpu, par("usr")[3] - 0.50, srt = 45, adj = 1, labels = xlab_cpu, xpd = NA)
    mtext("Non-NGS projects", side=3, line=-0.2, font=2)
    mtext(text="Mean CPU usage",side=1,line=2)

    my.plot.histogram(hd.oth_mem, xlim=xlim.mem, ylim=ylim.mem, col=oth.col, border=oth.col, lwd=hist.lwd, lend=hist.lend, xlab="", ylab="", axes=FALSE, main="")
    if (y.lines) segments(xlim.mem[1]-1, yat, xlim.mem[2], yat, lwd=0.5, col='grey', lty='dotted')
    axis(2, at=yat, labels=rep("", length(yat))) # y axis
    axis(1, at=xat_mem, labels=FALSE, mgp=c(1.5, 0.4, 0.1))
    #mtext(bquote(italic(N)[total]==.(n.oth)*','~italic(N)[scaled]==10^.(log10(scale.counts))), side=3, adj=1, cex=0.8, line=-0.2)
    mtext(bquote(italic(N)[total]==.(n.oth)~'(unscaled)'), side=3, adj=1, cex=0.8, line=-0.2)
    text(xat_mem, par("usr")[3] - 0.50, srt = 45, adj = 1, labels = xlab_mem, xpd = NA)
    mtext(text="Max GiB of RAM used",side=1,line=2)

    mtext(text=bquote('Job count in 2016, scaled sum'==10^.(log10(scale.counts))),side=2,line=0.6,outer=TRUE, las=0)

    if (do.pdf) dev.off()
}

my.plot.histogram <- function (x, freq = equidist, density = NULL, angle = 45, col = NULL, 
    border = par("fg"), lty = NULL, main = paste("Histogram of", 
        paste(x$xname, collapse = "\n")), sub = NULL, xlab = x$xname, 
    ylab, xlim = range(x$breaks), ylim = NULL, axes = TRUE, labels = FALSE, 
    add = FALSE, ann = TRUE, lwd=par("lwd"), lend=par("lend"), ...) 
{
    equidist <- if (is.logical(x$equidist)) 
        x$equidist
    else {
        h <- diff(x$breaks)
        diff(range(h)) < 1e-07 * mean(h)
    }
    if (freq && !equidist) 
        warning("the AREAS in the plot are wrong -- rather use 'freq = FALSE'")
    y <- if (freq) 
        x$counts
    else x$density
    nB <- length(x$breaks)
    if (is.null(y) || 0L == nB) 
        stop("'x' is wrongly structured")
    dev.hold()
    on.exit(dev.flush())
    if (!add) {
        if (is.null(ylim)) 
            ylim <- range(y, 0)
        if (missing(ylab)) 
            ylab <- if (!freq) 
                "Density"
            else "Frequency"
        plot.new()
        plot.window(xlim, ylim, "", ...)
        if (ann) 
            title(main = main, sub = sub, xlab = xlab, ylab = ylab, 
                ...)
        if (axes) {
            axis(1, ...)
            axis(2, ...)
        }
    }
    rect(x$breaks[-nB], 0, x$breaks[-1L], y, col = col, border = border, 
        angle = angle, density = density, lty = lty, lwd=lwd, lend=lend)
    if ((logl <- is.logical(labels) && labels) || is.character(labels)) 
        text(x$mids, y, labels = if (logl) {
            if (freq) 
                x$counts
            else round(x$density, 3)
        }
        else labels, adj = c(0.5, -0.5))
    invisible()
}

