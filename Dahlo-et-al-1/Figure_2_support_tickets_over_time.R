# Plot Figure 2 from Dahlö et al, containing support ticket counts, unique
# submitters, and tickets per submitter over time, using data from the
# tickets.dat file, "Figure_2_data_support_tickets.txt".

figure.type = 'pdf'
figure.create = switch(figure.type, pdf = grDevices::pdf,
                                    eps = grDevices::postscript,
                                    default = stop("unknown figure.type '", figure.type, "'", call.=FALSE))


study.end.date = "2017-01-01"
study.end.Date = as.Date(as.POSIXct(study.end.date, tz="GMT"))
study.min.duration = 60
output.file.tag = paste0(format(Sys.time(), '%Y%m%d'), '_', study.end.date)

options(prompt = "tix> ", stringsAsFactors = FALSE)
#reload = function(doit = FALSE) if (doit) source("Figure_2.R")
library(RColorBrewer)

# Figure sizes
# 85 mm half-page
half.width = 85 / 25.4
full.width = 170 / 25.4
short.height = 60 / 25.4
short50.height = 50 / 25.4
short40.height = 40 / 25.4
full.height = 75 / 25.4
tall.height = 120 / 25.4

# 170 mm full-page
# 225 max height figure and legend
# 300 dpi

#######
#######
####### support tickets
#######
#######
#
# Produce 3-panel plot of support ticket submission
#

tickets.dat = "Figure_2_data_support_tickets.txt"
plot.tickets.data = function(file = tickets.dat, do.pdf=FALSE) {
    dat = read.delim(file)
    dat$tu_ngs = dat$t_ngs / dat$u_ngs
    dat$tu_gen = dat$t_gen / dat$u_gen

    sink('fig_tickets.log', split=TRUE)

    cat('ticket dat: \n')
    print(dat)
    cat('sums: \n')
    print(colSums(dat[,2:7]))

    xlim = range(dat$year)
    ylim.a = with(dat, range(0, t_ngs,  t_gen)); ylim.a[2] = round(ylim.a[2]+500, -3)
    ylim.b = with(dat, range(0, u_ngs,  u_gen)); ylim.b[2] = round(ylim.b[2]+100, -3)
    ylim.c = with(dat, range(1, 5, tu_ngs, tu_gen))
    cat('ylim.a = ', ylim.a, '\n')
    cat('ylim.b = ', ylim.b, '\n')
    cat('ylim.c = ', ylim.c, '\n')

    if (do.pdf)
        figure.create(file=paste0('fig_tickets.',figure.type), width=half.width, height=short.height, pointsize=10)
    opa=par(mfcol=c(3, 1), mar=c(1.5, 4.6, 0.5, 0.5), oma=c(0, 0.2, 0.2, 0),
            lwd=2, mgp=c(1.9, 0.4, 0), tcl=-0.2, las=1, ps=10, cex.axis=0.8)

    ngs.col = brewer.pal(9, "Greens")[7]
    gen.col = brewer.pal(9, "Reds")[7]
    ngs.lty = 1
    gen.lty = 2
    y.ln.1 = 3.2
    y.ln.2 = 2.2
    y.ln.3 = 1.2
    y.cex = 0.80
    grid.lty = 3
    grid.col = brewer.pal(9, "Greys")[4]
    grid.lwd = 0.5
    do.grid = function(hh) {
        for (h in hh) segments(xlim[1]-0.5, h, xlim[2], h, lty=grid.lty, lwd=grid.lwd, col=grid.col)
    }
    do.tag = function(tag, offset) {
        x0 = par("usr")[1]; x1 = par("usr")[2]; y0 = par("usr")[3]; y1 = par("usr")[4]
        text(x = x0+((x1 - x0) * -0.165), y = y1 + ((y1 - y0) * +0.02), tag, font=2, cex=1.2, xpd=NA)
    }

    # panel A)  ticket numbers
    plot.new()
    plot.window(bty="n", xlim=xlim, ylim=ylim.a)
    do.grid(seq(ylim.a[1], ylim.a[2], by=1000))
    legend("topleft", inset=c(0.0, -0.15), legend=c("NGS projects", "Non-NGS projects"), 
           lty=c(ngs.lty, gen.lty), lwd=2, col=c(ngs.col, gen.col), y.intersp=0.9, cex=1.0, seg.len=2.4, xpd=NA,
           bty="o", box.col="white", bg="white")
    lines(dat$year, dat$t_gen, lty=gen.lty, lwd=2, col=gen.col)
    lines(dat$year, dat$t_ngs, lty=ngs.lty, lwd=2, col=ngs.col)
    axis(1, at=dat$year, tcl=-0.15, labels=FALSE)
    axis(2, at=seq(ylim.a[1], ylim.a[2], by=1000))
    mtext(side=2, line=y.ln.1, text="Tickets", las=0, cex=y.cex)
    mtext(side=2, line=y.ln.2, text="submitted", las=0, cex=y.cex)
    do.tag("A")

    # panel B)  ticket submitters
    plot.new()
    plot.window(bty="n", xlim=xlim, ylim=ylim.b)
    do.grid(seq(ylim.b[1], ylim.b[2], by=200))
    lines(dat$year, dat$u_gen, lty=gen.lty, lwd=2, col=gen.col)
    lines(dat$year, dat$u_ngs, lty=ngs.lty, lwd=2, col=ngs.col)
    axis(1, at=dat$year, tcl=-0.15, labels=FALSE)
    axis(2, at=seq(ylim.b[1], ylim.b[2], by=200))
    mtext(side=2, line=y.ln.1, text="Unique", las=0, cex=y.cex)
    mtext(side=2, line=y.ln.2, text="submitters", las=0, cex=y.cex)
    do.tag("B")

    # panel C)  tickets/submitter
    plot.new()
    plot.window(bty="n", xlim=xlim, ylim=ylim.c)
    do.grid(seq(ylim.c[1], ylim.c[2], by=1))
    lines(dat$year, dat$tu_gen, lty=gen.lty, lwd=2, col=gen.col)
    lines(dat$year, dat$tu_ngs, lty=ngs.lty, lwd=2, col=ngs.col)
    axis(1, at=dat$year, tcl=-0.15, labels=dat$year, cex.axis=1.1)
    axis(2, at=seq(ylim.c[1], ylim.c[2], by=1))
    #mtext(side=2, line=y.ln.1, text="Tickets", las=0, cex=y.cex)
    mtext(side=2, line=y.ln.2, text="Tickets per", las=0, cex=y.cex)
    mtext(side=2, line=y.ln.3, text="submitter", las=0, cex=y.cex)
    do.tag("C")

    par(opa)
    if (do.pdf) dev.off()

    sink()
    invisible(dat)
}
