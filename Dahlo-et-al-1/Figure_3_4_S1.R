# vim: set number:

# TODO: stats for unincluded jobs --

# 1. COMPLETED < 60 s
#
#     See totaltime.dat and the file totaltime_20151127_2015-11-01.txt (one example)
#
# 2. Same totaltime.dat for non-completed jobs
#
# 3. Plot/table of job class and duration statistics
#

fig.type = ifelse(exists('fig.type'), fig.type, 'pdf')
fig.create = switch(fig.type, pdf = grDevices::pdf,
                              eps = grDevices::postscript,
                              default = stop("unknown fig.type '", fig.type, "'", call.=FALSE))


study.end.date = "2017-01-01"
study.end.Date = as.Date(as.POSIXct(study.end.date, tz="GMT"))
study.min.duration = 60
output.file.tag = paste0(format(Sys.time(), '%Y%m%d'), '_', study.end.date)
source('legendx.R')


cat('

Producing figures in', fig.type, 'format

To generate figures and stats for all data:

   master.do.all(include.all.data = TRUE)

This does not include plot.install.data(do.pdf = TRUE)

')



# This is the workflow, but is for how-to. The R environment already contains
# jdat, jsdat, tab.ps, tab.bbs, jtabs.c, jtabs.ct

master.do.all = function(include.all.data = TRUE) {
    sink(paste0(output.file.tag, '.log'), split = TRUE)

    jdat = master.prepare.jobstate(all.data=TRUE, min.duration=1, end.date=study.end.date)
    jdat <<- jdat
    jsdat = jdat %>% filter(duration_sec >= 60)
    jsdat <<- jsdat
    prepare.jobstate.tables(jdat, min.duration = 60)  # creates tab.ps and tab.bbs
    tab.ps <<- tab.ps
    tab.bbs <<- tab.bbs
    plot.jobstate.trends(do.pdf=TRUE)  # uses the above tables
    jtabs.c = generate.projgroup.tables(jdat, include.timeout=FALSE)
    jtabs.ct = generate.projgroup.tables(jdat, include.timeout=TRUE)
    jtabs.c <<- jtabs.c
    jtabs.ct <<- jtabs.ct
    all.plot.fig.projgroup(do.pdf=TRUE, include.timeout=TRUE)
    sink()
    "master.do.all completed"
}

#--------------------------------


options(prompt = "Figure_3_4_S1> ", stringsAsFactors = FALSE)
reload = function(doit = FALSE) if (doit) source("Figure_3_4_S1")
library(dplyr)
library(stringr)
library(lattice)
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
####### Module and library installs
#######
#######
#
# Produce plot of module installs and updates for bio and non-bio projects
#
# categories are b, o, s
# b is bioinformatics software
# o is nonbioinformatics specialty software or software that extends capabilities such as general-purpose libraries, solvers or debuggers
# s is system software including compilers and interpreters, and other non-b, non-o events including the installation of data sources
#
# note we installed the 1000-genomes datasets, Simons, Blast databases, etc.

install.dat.file = "Figure_3_data.txt"
install.type = c("oi", "ou", "bi", "bu")
install.type.long = setNames(c("Non-NGS new installs", "Non-NGS updates",
                               "NGS-related new installs", "NGS-related updates"),
                             install.type)
install.type.cols = brewer.pal(8, "Paired")[c(7,8,3,4)]; names(install.type.cols) = install.type
plot.install.data = function(file = install.dat.file, end.date = study.end.date, start.date = "2014-01-01", do.pdf=FALSE) {
    dat = tbl_df(read.delim(file, sep="\t", colClasses=c("character","character","character","character"),
                            header=TRUE, fill=TRUE))
  
    dat = dat %>%
        mutate(date = as.Date(as.POSIXct(date, tz="GMT"))) %>%
        mutate(install_type = paste0(cat, action)) %>%
        select(-cat, -action, -module) %>%
        filter(install_type %in% install.type) %>%
        block.data(blocks.1month.2014.2016, "date", end.date = end.date, start.date = start.date)
    start.date = as.Date(as.POSIXct(start.date, tz="GMT"))
    end.date = as.Date(as.POSIXct(end.date, tz="GMT"))
    dat = dat %>% 
        filter(date >= start.date & date < end.date)
    d = with(dat, table(date_block, install_type))
    d = apply(d, 2, cumsum)
    d = d.orig = d[, install.type]
    dplot = data.frame(date_block = rownames(d), t(apply(d, 1, cumsum)))
    # 
    xpos = setNames(seq(along=blocks.1month.2014.2016), blocks.1month.2014.2016)
    xpos = xpos[-length(xpos)]
    dplot$xpos = xpos
    this.ylim = range(dplot[, 2:5])
    this.ylim = range(0, ceiling(this.ylim[2]/100)*100)
    if (do.pdf)
        fig.create(file=paste0("fig_installs.", fig.type), width=half.width, height=short50.height, pointsize=10)
    opa=par(mar=c(2.8, 2.9, 0.7, 0), lwd=2, mgp=c(1.9, 0.4, 0), yaxs='i', tcl=-0.2, xpd=NA, ps=9, cex.axis=0.9)
    plot.new()
    plot.window(bty="n", xlim=range(xpos), ylim=this.ylim)
    # grid lines
    yat = axTicks(2)
    abline( h = yat, lty = 3, lwd=1.0, col = brewer.pal(9, "Greys")[4], xpd=FALSE)
    # axes
    axis(1, at=xpos, tcl=-0.15, labels=FALSE)
    blab = blocks.to.label[blocks.to.label %in% blocks.1month.2014.2016]
    tck = axis(1, at=xpos[blab], tcl=-0.3, labels=FALSE)
    text(tck, 1, labels=blocks.shortlabels[blab], srt=45, adj=c(1.1, 1.65), cex=0.9)
    axis(2, at=yat, las=2, cex=0.8)
    # installations and updates
    for (it in 2:5) {
        lines(dplot$xpos, dplot[, it], col=install.type.cols[names(dplot)[it]], lty=1)
    }
    mtext("Research software installs & updates", side=2, line=1.9, cex=0.9, adj=0.8)
    lgnd = rev(paste0(install.type.long, " (", d.orig[nrow(d.orig),], ")"))
    # legendx() !!!!
    legendx("topleft", inset=c(0.008, -0.15), cex=0.9,
            fill=rev(install.type.cols), border=rev(install.type.cols),
            y.intersp=0.8,
            bty="o", box.col='white', bg='white', box.hex=0.8,
            legend=lgnd)
    par(opa)
    if (do.pdf) dev.off()
    dplot$n_oi = dplot$oi
    dplot$n_ou = dplot$ou - dplot$oi
    dplot$n_bi = dplot$bi - dplot$ou
    dplot$n_bu = dplot$bu - dplot$bi
    sink("fig_installs.txt")
    print(dplot)
    sink()
    dplot
}


#---- clusters
default.clusters = c("milou", "kalkyl", "tintin"); names(default.clusters) = default.clusters
valid.clusters = c(default.clusters, "halvan", "fysast1", "nestor"); names(valid.clusters) = valid.clusters
cluster.cols = brewer.pal(length(valid.clusters), "Set1"); names(cluster.cols) = valid.clusters
cluster.colours = function(pr) {
  stopifnot(all(pr %in% valid.clusters))
  cluster.cols[pr]
}
get.cores = function(cl=default.clusters) {
  switch(match.arg(cl, valid.clusters), milou = 16, kalkyl = 8, tintin = 16, halvan = 64, fysast1 = 16, nestor = 16)
}

#---- projects
valid.projgroups = c("p+s", "b+bs", "a", "g")
projgroup.type = function(pr=valid.projgroups) {
  switch(match.arg(pr),
         'b+bs' = "NGS", 
         'p+s'  = "Non-NGS", 
         'a'    = "Sequencing platform", 
         'g'    = "Course", 
         default="Unknown")
}
projgroup.name = function(pr=valid.projgroups) {
  ans = projgroup.type(pr); if(pr != "a") ans = paste(ans, "projects"); ans
}
projgroup.clust = list('b+bs' = c("kalkyl","milou","tintin"), 
                       'p+s'  = c("kalkyl","tintin","fysast1"), 
                       'a'    = c("kalkyl","milou","nestor"), 
                       'g'    = c("kalkyl","milou","tintin")) 
stopifnot(all(sort(names(projgroup.clust)) == sort(valid.projgroups)))  # make sure all projects accounted for
projgroup.numclust = unlist(lapply(projgroup.clust, length))
projgroup.clusters = function(pr=valid.projgroups) {
  projgroup.clust[[match.arg(pr)]]
}
projgroup.legend.clusters = function(pr=valid.projgroups) {
    legend.clust = list('b+bs' = c("kalkyl","tintin","fysast1","milou"), 
                        'p+s'  = NULL,
                        'a'    = c("kalkyl","milou","nestor"), 
                        'g'    = c("kalkyl","milou","tintin")) 
    legend.clust[[match.arg(pr)]]
}

valid.projects = c("b","a","p","s","g","bs")
project.type = function(pr=valid.projects) {
  switch(match.arg(pr),
         a = "Sequencing platform", 
         b = "NGS", 
         g = "Course", 
         p = "Non-NGS", 
         s = "SNIC Large-Scale", 
         bs = "SNIC NGS", 
         default="Unknown")
}
project.name = function(pr=valid.projects) {
  ans = project.type(pr); if(pr != "a") ans = paste(ans, "projects"); ans
}
project.clust = list(b=c("kalkyl","milou"), 
                     a=c("kalkyl","milou","nestor"), 
                     p=c("kalkyl","tintin","fysast1"), 
                     s=c("kalkyl","tintin"), 
                     g=c("kalkyl","milou","tintin"), 
                     bs=c("kalkyl","tintin")) 
stopifnot(all(sort(names(project.clust)) == sort(valid.projects)))  # make sure all projects accounted for
project.numclust = unlist(lapply(project.clust, length))
project.clusters = function(pr=valid.projects) {
  project.clust[[match.arg(pr)]]
}

#---- job types
valid.jobtype = c("core", "partial", "node", "multi")
ltys.jobtype = c("dotted", "dashed", "solid", "81"); names(ltys.jobtype) = valid.jobtype
jobtype.ltys = function(jt) {
  stopifnot(all(jt %in% valid.jobtype))
  ltys.jobtype[jt]
}
cluster.jobtype.colours = function(cl=valid.clusters) {
  cols = switch(match.arg(cl),
                kalkyl  = brewer.pal(9, "Blues")  [c(2,4,6,8)], #6,7)],
                milou   = brewer.pal(9, "Greens") [c(2,4,6,8)], #6,7)],
                tintin  = brewer.pal(9, "Reds")   [c(2,4,6,8)], #6,7)],
                fysast1 = brewer.pal(9, "Purples")[c(2,4,6,8)], #6,7)],
                nestor  = brewer.pal(9, "Purples")[c(2,4,6,8)], #6,7)],
                halvan  = brewer.pal(9, "Oranges")[c(2,4,6,8)]) #6,7)])
  names(cols) = valid.jobtype
  cols
}


#---- blocks and labels for blocks
blocks.4month = c("2010-10-01", 
                  "2011-02-01", "2011-06-01", "2011-10-01", 
                  "2012-02-01", "2012-06-01", "2012-10-01", 
                  "2013-02-01", "2013-06-01", "2013-11-15", 
                  "2014-03-01", "2014-07-01", "2014-10-01",
                  "2015-02-01", "2015-06-01", "2015-10-01",
                  "2016-02-01", "2016-06-01", "2016-10-01",
                  "2016-02-01") # last entry should always be in the future
dates.4month = as.Date(as.POSIXct(blocks.4month, tz="GMT"))
names(dates.4month) = names(blocks.4month) = blocks.4month
last.date.4month = dates.4month[length(dates.4month)]

blocks.1month.2014.2016 = c(
                  "2014-01-01", "2014-02-01", "2014-03-01", "2014-04-01", "2014-05-01", "2014-06-01",
                  "2014-07-01", "2014-08-01", "2014-09-01", "2014-10-01", "2014-11-01", "2014-12-01",
                  "2015-01-01", "2015-02-01", "2015-03-01", "2015-04-01", "2015-05-01", "2015-06-01",
                  "2015-07-01", "2015-08-01", "2015-09-01", "2015-10-01", "2015-11-01", "2015-12-01",
                  "2016-01-01", "2016-02-01", "2016-03-01", "2016-04-01", "2016-05-01", "2016-06-01",
                  "2016-07-01", "2016-08-01", "2016-09-01", "2016-10-01", "2016-11-01", "2016-12-01",
                  "2017-01-01") # last entry should always be in the future
blocks.1month = c("2010-10-01", "2010-11-01", "2010-12-01",
                  "2011-01-01", "2011-02-01", "2011-03-01", "2011-04-01", "2011-05-01", "2011-06-01",
                  "2011-07-01", "2011-08-01", "2011-09-01", "2011-10-01", "2011-11-01", "2011-12-01",
                  "2012-01-01", "2012-02-01", "2012-03-01", "2012-04-01", "2012-05-01", "2012-06-01",
                  "2012-07-01", "2012-08-01", "2012-09-01", "2012-10-01", "2012-11-01", "2012-12-01",
                  "2013-01-01", "2013-02-01", "2013-03-01", "2013-04-01", "2013-05-01", "2013-06-01",
                  "2013-07-01", "2013-08-01", "2013-09-01", "2013-10-01", "2013-11-01", "2013-12-01",
                  "2014-01-01", "2014-02-01", "2014-03-01", "2014-04-01", "2014-05-01", "2014-06-01",
                  "2014-07-01", "2014-08-01", "2014-09-01", "2014-10-01", "2014-11-01", "2014-12-01",
                  "2015-01-01", "2015-02-01", "2015-03-01", "2015-04-01", "2015-05-01", "2015-06-01",
                  "2015-07-01", "2015-08-01", "2015-09-01", "2015-10-01", "2015-11-01", "2015-12-01",
                  "2016-01-01", "2016-02-01", "2016-03-01", "2016-04-01", "2016-05-01", "2016-06-01",
                  "2016-07-01", "2016-08-01", "2016-09-01", "2016-10-01", "2016-11-01", "2016-12-01",
                  "2017-01-01") # last entry should always be in the future
dates.1month = as.Date(as.POSIXct(blocks.1month, tz="GMT"))
names(dates.1month) = names(blocks.1month) = blocks.1month
last.date.1month = dates.1month[length(dates.1month)]

blocks.to.label = c("2009-12-01", 
                    "2010-12-01", "2011-03-01", "2011-06-01", "2011-09-01",
                    "2011-12-01", "2012-03-01", "2012-06-01", "2012-09-01",
                    "2012-12-01", "2013-03-01", "2013-06-01", "2013-09-01",
                    "2013-12-01", "2014-03-01", "2014-06-01", "2014-09-01",
                    "2014-12-01", "2015-03-01", "2015-06-01", "2015-09-01",
                    "2015-12-01", "2016-03-01", "2016-06-01", "2016-09-01",
                    "2016-12-01")
blocks.shortlabels = substr(blocks.to.label, 1, 7); names(blocks.shortlabels) = blocks.to.label

# collective globals for the plotting functions
# drop the unused blocks from the set to use
used.blocks = blocks.1month[dates.1month < study.end.Date]; names(used.blocks) = used.blocks
blocks.to.label = blocks.to.label[blocks.to.label %in% used.blocks]  # shrink set down to those used
all.block.labels = used.blocks; all.block.labels[! all.block.labels %in% blocks.to.label] = ""
all.block.shortlabels = substr(all.block.labels, 1, 7)
xpos = seq_along(used.blocks); names(xpos) = used.blocks

#-------------------------- summarise data --------------


generate.tables = function(dat) {
    ans = list()
    attach(dat, warn.conflicts=FALSE)
    ans[["jobs-by-cluster"]] = table(cluster)
    ans[["jobs-by-cluster,proj"]] = table(cluster, proj)
    ans[["jobs-by-cluster:proj,year"]] = table(paste(sep=":", cluster, proj), format(start_date, "%Y"))
    ans[["jobs-by-cluster:proj,date_block"]] = table(paste(sep=":", cluster, proj), date_block)
    ans[["jobs-by-cluster:proj:jobtype,date_block"]] = table(paste(sep=":", cluster, proj, job_type), date_block)
    ans[["df-core_hours-cluster,proj,job_type,date_block"]] =  generate.dataframe(dat, c("cluster","proj","job_type","date_block"), "core_hours")
    ans[["df-job_count-cluster,proj,job_type,date_block"]] =  generate.dataframe(dat, c("cluster","proj","job_type","date_block"), "job_count")
    sapply(names(ans), function(x) { cat("\n\n",x,"\n"); print(ans[[x]]) })
    detach(dat)
    invisible(ans)
}

all.plot.fig.projgroup = function(do.pdf=TRUE, include.timeout=TRUE) {
    all.plot.fig.projgroup.grp(1, do.pdf, include.timeout)
    all.plot.fig.projgroup.grp(2, do.pdf, include.timeout)
    plot.fig.corehours.projgroup(pr='a', all.data=TRUE, do.pdf=do.pdf)  #includes timeout
}

all.plot.fig.projgroup.grp = function(grp, do.pdf=TRUE, include.timeout=FALSE) {
    if (include.timeout) {
        dat = jtabs.ct[["df-core_hours-cluster,projgroup,job_type,date_block"]]
        jobstate = "Completed and Timeout"
        filetag = paste0(grp, "_ct")
    } else {
        dat = jtabs.c[["df-core_hours-cluster,projgroup,job_type,date_block"]]
        jobstate = "Completed"
        filetag = paste0(grp, "_c")
    }
    if (do.pdf)
        fig.create(file=paste0("fig_corehours-projgroup_",filetag,".",fig.type), width=full.width, height=90/25.4)
    opa=par(mfcol=c(2, 1), oma=c(0.8, 0.3, 0, 0), mar=c(1.0, 1.5, 1.1, 0), lwd=2, mgp=c(1.0, 0.3, 0), tcl=-0.2, xpd=NA, ps=9, cex.axis=0.7) #0.8
    plot.fig.corehours.projgroup(dat, tag="A) ", pr=ifelse(grp==1,'b+bs', 'a'),
                                 all.data=TRUE, do.pdf=FALSE, add=TRUE, js=jobstate)
    tck = plot.fig.corehours.projgroup(dat, tag="B) ", pr=ifelse(grp==1, 'p+s', 'g'),
                                       all.data=TRUE, do.pdf=FALSE, add=TRUE, js=jobstate)
    text(tck, 1, labels=all.block.shortlabels, srt=45, adj=c(1.1, 1.2), cex=0.7)
    mtext("Monthly core hours booked (millions)", side=2, outer=TRUE, line=-0.7, cex=0.9, las=0)
    if (do.pdf) dev.off()
    par(opa)
}



# by default, Completed and Timeout jobs
plot.fig.corehours.projgroup = function(dat=jtabs.ct[["df-core_hours-cluster,projgroup,job_type,date_block"]],
                              pr=valid.projgroups,
                              cl=c("milou"),
                              jt=valid.jobtype,
                              all.data=all.data,
                              tag="",
                              do.pdf=FALSE,
                              add=FALSE,
                              js="Completed and Timeout",
                              y.lines=TRUE) {
    if (add)               { ht = 40/25.4; s1 = 0.5; xal = FALSE }
    else if (tag == "A) ") { ht = 35/25.4; s1 = 0.5; xal = FALSE }
    else if (tag == "B) ") { ht = 45/25.4; s1 = 1.5; xal = TRUE }
    else                   { ht = 60/25.4; s1 = 1.5; xal = TRUE }
    pr = match.arg(pr)
    if (missing(cl)) cl = projgroup.clusters(pr)
    stopifnot(all(jt %in% valid.jobtype))
    if (pr %in% c('a','g'))
         this.ylim = c(0, 600000)
    else this.ylim = c(0, 2500000)
    dat = dat %>%
        mutate(cl.jt = paste(sep=":", cluster, job_type)) %>%
        filter(cluster %in% cl, projgroup %in% pr) %>%
        filter(as.Date(as.POSIXct(date_block, tz="GMT")) < last.date.1month) %>%  # remove last month
        filter(! (cluster == "kalkyl" & date_block == "2014-01-01"))     # remove last month of kalkyl operation
    tab = xtabs(core_hours ~ cl.jt + date_block, dat) 
    cl.jt.ord = paste(sep=":", rep(cl, rep(length(jt), length(cl))), rep(jt, length(cl)))
    # tab might be missing rows from cl.jt.ord, add them if so
    m = cl.jt.ord[! cl.jt.ord %in% rownames(tab)]
    z = matrix(0, length(m), ncol(tab), dimnames=list(column=m, date_block=dimnames(tab)[[2]]))
    tab = rbind(tab, z)
    # now reorder tab
    tab = tab[cl.jt.ord, ]  # cluster order, then job booking order: core < partial < node < multi
    tab.cols = unlist(lapply(cl, cluster.jobtype.colours))
    if (do.pdf)
        fig.create(file=paste0("fig_corehours-projgroup_", pr, "_", paste(collapse="+", cl),
                        ifelse(all.data, "_alldata", "_subset"), ".", fig.type), width=full.width, height=ht)
    if (! add)
        opa=par(mar=c(s1, 1.5, 1.6, 0), lwd=2, mgp=c(1.0, 0.3, 0), tcl=-0.2, xpd=NA, ps=9, cex.axis=0.7)
    if (add && (tag == 'B) ' || tag == 'B) '))
        tck = barplot(tab, space=-0.03, ylim=this.ylim, border=NA, col=tab.cols, axes=FALSE, axisnames=FALSE, las=3)
    else tck = barplot(tab, space=-0.03, ylim=this.ylim, border=NA, col=tab.cols, axes=FALSE, axisnames=FALSE, las=3)
    usr = par('usr')
    if (y.lines) {
        y.lines = pretty(this.ylim)
        segments(-2, y.lines, ncol(tab)-2, y.lines, lwd=0.5, col="grey", lty='dotted')
    }
    #axis(1, at=tck, labels=all.block.shortlabels, las=3)
    if (xal)
        text(tck, 1, labels=all.block.shortlabels, srt=45, adj=c(1.0, 1.0), cex=0.7)
    else if (add) {        # && tag == "A) ") {
        at = tck[all.block.shortlabels != ""]
        sub.at = tck[all.block.shortlabels == ""]
        #axis(1, at = at, labels=NA)
        l = (usr[4] - usr[3]) * 0.02  # ticks will be 2% of plot y extent
        segments(at, usr[3], at, usr[3]-l, lwd=1.0)
        l = (usr[4] - usr[3]) * 0.01  # subticks will be 1% of plot y extent
        segments(sub.at, usr[3], sub.at, usr[3]-l, lwd=0.5)
    }
    segments(0, usr[3], ncol(tab)-2, usr[3], lwd=0.5, col="black", lty=1)
    yat = axTicks(2)
    rnd = function(x) {
        if (!x) "0"
        else if (max(yat) > 500000) sprintf("%.1f", x)
        else if (max(yat) > 100000) sprintf("%.1f", x)
        else sprintf("%.2f", x)
    }
    axis(2, las=2, line=-0.5, at=yat, labels=sapply(yat/(10^6), rnd), cex=0.7)
    if (! add)
        mtext('Monthly core hours booked (millions)', side=2, line=0.7, cex=0.9)
    if (add && tag == "B) ") {
        # do not produce a legend
        lgnd = NULL
    } else if (add && tag == "A) ") {
        lcl = projgroup.legend.clusters(pr)
        lgnd = paste(sep=" ", rep(lcl, rep(length(jt), length(lcl))), rep(jt, length(lcl)))
        tab.cols = unlist(lapply(lcl, cluster.jobtype.colours))
    } else {
        lgnd = str_replace(cl.jt.ord, ":", " ")
        lcl = cl
    }
    #cat("strwidth = ", strwidth(lgnd), "\n")
    #cat("max strwidth = ", max(strwidth(lgnd)), "\n")
    if (! is.null(lgnd))
        # NOTE: legendx() !!!
        legendx("topleft", inset=c(0.02, 0.00),
               title="Cluster and job type", title.adj=0.05, title.font=2,
               ncol=length(lcl),
               bty="o", box.col="white", bg="white",
               col=tab.cols, fill=tab.cols, border=tab.cols,
               cex=0.7, x.intersp=0.2, y.intersp=0.8, 
               legend=lgnd, text.width=5.5)
    if (add) mtext(paste0(tag, projgroup.name(pr), ": ", js, " jobs"), font=2, side=3, line=+0.3,
                   adj=-0.02, las=0, cex=1.0)
    else mtext(paste0(tag, projgroup.name(pr), ": ", js, " jobs"), font=2, side=3, line=-0.8,
               outer=TRUE, adj=0.05, las=0, cex=1.0)
    #mtext("Completed jobs", font=2, side=3, line=-1.6, outer=TRUE, adj=0.02, las=0, cex=1.0)
    if (! add) par(opa)
    if (do.pdf) dev.off()
    ####
    if (add) tck
    else invisible(list(projgroup=pr, clusters=cl, jobtypes=jt, table=tab, colours=tab.cols))
}

generate.projgroup.tables = function(dat, include.timeout=FALSE) {
    ans = list()
    if (include.timeout) 
        dat = filter(dat, jobstate %in% c('COMPLETED', 'TIMEOUT'), duration_sec >= 60)
    else dat = filter(dat, jobstate == 'COMPLETED', duration_sec >= 60)
    attach(dat, warn.conflicts=FALSE)
    ans[["df-core_hours-cluster,projgroup,job_type,date_block"]] =  generate.dataframe(dat, c("cluster","projgroup","job_type","date_block"), "core_hours")
    ans[["df-job_count-cluster,projgroup,job_type,date_block"]] =  generate.dataframe(dat, c("cluster","projgroup","job_type","date_block"), "job_count")
    sapply(names(ans), function(x) { cat("\n\n",x,"\n"); print(ans[[x]]) })
    detach(dat)
    invisible(ans)
}

generate.dataframe = function(dat, type, response=c("job_count","core_hours")) {
    response = match.arg(response)
    dat = tbl_df(dat) %>% 
        mutate(year = format(start_date, "%Y")) %>%
        mutate(core_hours = cores * duration_sec / 3600.) %>%
        mutate(date_block = ordered(date_block)) %>%
        group_by_(.dots = type)
    stopifnot(all(type %in% names(dat)))
    if (response == "job_count") {
        dat = dat %>% summarise(job_count = length(cluster))
    } else if (response == "core_hours") {
        dat = dat %>% summarise(core_hours = sum(core_hours))
    } else {
        stop(response,"is an unknown response in generate.dataframe")
    }
    dat
}


#-------------------------- prepare data --------------

# # The steps required to download a new job database and prepare the data in R
#
# # In the shell, in the data/ directory:
# ./update-dbs.sh    # this requires entry of my uppmax password
# ./extract-from-db.sh  # this version implements separate extraction of b-type s projects
#
# # Then, within R in the base directory:
# reload(T)  # just to make sure you have the latest of everything
# dat = master.prepare.data()
# tabs = generate.tables(dat)
# # Generate the figures
# corehours.stats = master.plot.corehours()
# jobcount.stats = master.plot.jobcount()

# To generate data and plots for *all* jobs, not just a 10% sample
# all.dat = master.prepare.data(TRUE)
# tabs = generate.tables(all.dat)  # this must be named 'tabs' for figure generation to work
# # Generate the figures
# corehours.stats = master.plot.corehours()
# jobcount.stats = master.plot.jobcount()


master.prepare.jobstate = function(all.data=FALSE, min.duration=1, end.date=study.end.date) {
    logfile = paste0("master.prepare.jobstate_all.data", all.data, "_end.date", end.date, ".log")
    sink(logfile, split = TRUE)
    cat("master.prepare.jobstate: writing log to", logfile, "\n")
    cat("master.prepare.jobstate: all.data =", all.data, "end.date =", end.date, "\n")
    subs = ifelse(all.data, 100, 10)
    do.prepare.data = function(cls, pr, sb, md = min.duration, subext="-COMPLETED") {
        prepare.data(cls, subsample.percent=sb,  min.duration=md, proj=pr, end.date=end.date, subextension=subext)
    }
    d1c  = do.prepare.data(c("kalkyl","milou"),  "b",  subs, subext="-COMPLETED")
    d1n  = do.prepare.data(c("kalkyl","milou"),  "b",  subs, subext="-NOTCOMPLETED")
    d2c  = do.prepare.data(c("kalkyl","milou","nestor"),  "a",  subs, subext="-COMPLETED")
    d2n  = do.prepare.data(c("kalkyl","milou","nestor"),  "a",  subs, subext="-NOTCOMPLETED")
    d3c  = do.prepare.data(c("tintin"),          "p",  subs, subext="-COMPLETED")
    d3n  = do.prepare.data(c("tintin"),          "p",  subs, subext="-NOTCOMPLETED")
    d4c  = do.prepare.data(c("tintin"),          "s",  subs, subext="-COMPLETED")
    d4n  = do.prepare.data(c("tintin"),          "s",  subs, subext="-NOTCOMPLETED")
    d5c  = do.prepare.data(c("fysast1"),         "p",  100,  subext="-COMPLETED")
    d5n  = do.prepare.data(c("fysast1"),         "p",  100,  subext="-NOTCOMPLETED")
    d6c  = do.prepare.data(c("halvan"),          "b",  100,  subext="-COMPLETED")
    d6n  = do.prepare.data(c("halvan"),          "b",  100,  subext="-NOTCOMPLETED")
    d7c  = do.prepare.data(c("halvan"),          "p",  100,  subext="-COMPLETED")
    d7n  = do.prepare.data(c("halvan"),          "p",  100,  subext="-NOTCOMPLETED")
    d8c  = do.prepare.data(c("kalkyl","milou","tintin","fysast1"), "g", 100, subext="-COMPLETED")
    d8n  = do.prepare.data(c("kalkyl","milou","tintin","fysast1"), "g", 100, subext="-NOTCOMPLETED")
    d9c  = do.prepare.data(c("halvan"),          "bs", 100,  subext="-COMPLETED")
    d9n  = do.prepare.data(c("halvan"),          "bs", 100,  subext="-NOTCOMPLETED")
    d10c = do.prepare.data(c("kalkyl","tintin"), "bs", 100,  subext="-COMPLETED")
    d10n = do.prepare.data(c("kalkyl","tintin"), "bs", 100,  subext="-NOTCOMPLETED")
    d11c = do.prepare.data(c("kalkyl"),          "p",  subs, subext="-COMPLETED")
    d11n = do.prepare.data(c("kalkyl"),          "p",  subs, subext="-NOTCOMPLETED")
    d12c = do.prepare.data(c("kalkyl"),          "s",  subs, subext="-COMPLETED")
    d12n = do.prepare.data(c("kalkyl"),          "s",  subs, subext="-NOTCOMPLETED")
    sink()
    #totaltime.dat <- rbind(f(d1), f(d2), f(d3), f(d4), f(d5), f(d6), f(d7), f(d8), f(d9), f(d10))
    p.s  = bind_rows(d3c, d3n, d4c, d4n, d5c, d5n, d7c, d7n, d11c, d11n, d12c, d12n);  p.s$projgroup = "p+s"
    b.bs = bind_rows(d1c, d1n, d6c, d6n, d9c, d9n, d10c, d10n); b.bs$projgroup = "b+bs"
    g    = bind_rows(d8c, d8n)                                ;    g$projgroup = "g"
    a    = bind_rows(d2c, d2n)                                ;    a$projgroup = "a"
    f = function(x) attr(x, "xdat")
    ttjd <- rbind(      f(d1c), f(d1n), f(d2c), f(d2n), f(d3c), f(d3n), f(d4c), f(d4n))
    ttjd <- rbind(ttjd, f(d5c), f(d5n), f(d6c), f(d6n), f(d7c), f(d7n), f(d8c), f(d8n))
    ttjd <- rbind(ttjd, f(d9c), f(d9n), f(d10c), f(d10n))
    ttjd <- rbind(ttjd, f(d11c), f(d11n), f(d12c), f(d12n))
    rownames(ttjd) = NULL
    write.data(ttjd, paste0('totaltime_jobstate_', output.file.tag, '.txt'), gzip = FALSE)
    totaltime.jobstate.dat <<- ttjd
    #bind_rows(d1, d2, d3, d4, d5, d6, d7, d8, d9, d10) %>% arrange(start_date)
    bind_rows(b.bs, p.s, g, a) %>% arrange(start_date)
}

prepare.jobstate.tables = function(dat, min.duration = 60) {
    dat = dat %>%
        select(proj_id, start_date, date_block, cluster, jobstate, job_type, duration_sec, core_hours, proj, projgroup)
    dat = dat %>%
        filter(duration_sec >= min.duration) %>%
        filter(jobstate != 'NODE_FAIL') %>%
        mutate(jobstate = ordered(jobstate, levels=c('COMPLETED','TIMEOUT','CANCELLED','FAILED')))
    dat.ps  = dat %>% filter(projgroup == 'p+s')
    dat.bbs = dat %>% filter(projgroup == 'b+bs')
    # job numbers
    # tab.n.ps  = with(dat.ps,  table(jobstate, date_block))
    # tab.n.bbs = with(dat.bbs, table(jobstate, date_block))
    # core hours
    tab.ps  <<- xtabs(core_hours ~ jobstate + date_block, dat.ps)
    tab.bbs <<- xtabs(core_hours ~ jobstate + date_block, dat.bbs)
    "tab.ps for p+s, tab.bbs for b+bs"
}

jobstate = c('COMPLETED', 'TIMEOUT', 'CANCELLED', 'FAILED', 'NODE_FAIL')
names(jobstate) = jobstate
jobstate.display.l = setNames(c("Completed", "Timeout", "Cancelled", "Failed", "Other"), jobstate)
jobstate.display = jobstate.display.l[-length(jobstate.display.l)] # drop NODE_FAIL/Other

#plot.jobstate.trends = function(tab.ps, tab.bbs, do.pdf=FALSE) {
plot.jobstate.trends = function(do.pdf=TRUE, y.lines=TRUE) {
    cat("tab.ps  tab.bbs  taken from the environment\n")
    #this.ylim = range(0, 2500000, apply(tab.ps, 2, sum), apply(tab.bbs, 2, sum))
    this.ylim = range(0, apply(tab.ps, 2, sum), apply(tab.bbs, 2, sum))
    cat("this.ylim =", this.ylim, "\n")
    col.ps  = brewer.pal(9, "YlOrBr")[c(3,5,7,9)]
    col.bbs = brewer.pal(9, "BuGn")  [c(3,5,7,9)]
    if (do.pdf)
        fig.create(file=paste0("fig_jobstate.", fig.type), width=full.width, height=short50.height)
    opa = par(mfcol=c(1, 2), oma=c(0, 1, 0.5, 0), mar=c(1.7, 0.5, 0.5, 0), lwd=2, 
              mgp=c(1.0, 0.3, 0), tcl=-0.2, xpd=NA, ps=9, cex.axis=0.7)
    # panel A: p+s
    do.panel = function(tag = "A)", this.tab, this.cols) {
        tck = barplot(this.tab, space=-0.03, border=NA, col=this.cols, axes=FALSE,
                      axisnames=FALSE, las=3, ylim=this.ylim)
        if (y.lines) {
            y.lines = pretty(this.ylim)
            segments(0, y.lines, ncol(this.tab)-2, y.lines, lwd=0.5, col='grey', lty='dotted')
        }
        this.x.labels = all.block.shortlabels
        this.x.labels[seq(2, length(this.x.labels), by=2)] = ""
        at = tck[all.block.shortlabels != ""]  # every third month
        sub.at = tck[all.block.shortlabels == ""]
        usr = par('usr')
        l = (usr[4] - usr[3]) * 0.02
        segments(at, usr[3], at, usr[3]-l, lwd=0.75)
        l = (usr[4] - usr[3]) * 0.01
        segments(sub.at, usr[3], sub.at, usr[3]-l, lwd=0.5)
        text(tck, 1, labels=this.x.labels, srt=45, adj=c(1.1, 1.1), cex=0.7)
        yat = axTicks(2)
        rnd = function(x) {
            if (!x) "0"
            else if (max(yat) > 500000) sprintf("%.1f", x)
            else if (max(yat) > 100000) sprintf("%.2f", x)
            else sprintf("%.3f", x)
        }
        axis(2, las=2, line=-0.5, at=yat, labels=sapply(yat/(10^6), rnd), cex=0.7)
        legendx("topleft", inset=c(0.025, -0.04), ncol=1,
                bty="n", box.lwd=0, box.col=NA, bg='white',
                col=this.cols, fill=this.cols, border=this.cols,
                cex=0.7, x.intersp=0.2, y.intersp=0.8, 
                legend=c(jobstate.display))
        segments(0, usr[3], ncol(this.tab)-2, usr[3], lwd=0.5, col="black", lty=1)
        mtext(tag, font=2, side=3, line=+0.2, adj=-0.06, las=0, cex=1.0)
    }
    do.panel("A) NGS projects: all jobs", tab.bbs, col.bbs)
    mtext("Monthly core hours booked (millions)", side=2, line=0.5, cex=0.9)
    # panel B: p+s
    do.panel("B) Non-NGS projects: all jobs", tab.ps, col.ps)
    # annotation
    par(opa)
    if (do.pdf) dev.off()
}

#    job_id  proj_id      start partition cores cluster nodes  jobstate start_date is_core is_node
#    (int)    (chr)      (int)     (chr) (dbl)   (chr) (dbl)     (chr)     (date)   (lgl)   (lgl)
#1  157687 b2010010 1285917136      core     1  kalkyl     1 COMPLETED 2010-10-01    TRUE   FALSE
#2  157714 b2010010 1285917136      core     1  kalkyl     1 COMPLETED 2010-10-01    TRUE   FALSE
#3  157802 b2010007 1285930774      node     8  kalkyl     1 COMPLETED 2010-10-01   FALSE    TRUE
#4  157820 b2010010 1285939193      core     1  kalkyl     1 COMPLETED 2010-10-01    TRUE   FALSE
#5  157824 b2010010 1285939193      core     1  kalkyl     1 COMPLETED 2010-10-01    TRUE   FALSE
#6  157825 b2010010 1285939205      core     1  kalkyl     1 COMPLETED 2010-10-01    TRUE   FALSE
#7  157827 b2010010 1285939205      core     1  kalkyl     1 COMPLETED 2010-10-01    TRUE   FALSE
#8  157828 b2010010 1285939205      core     1  kalkyl     1 COMPLETED 2010-10-01    TRUE   FALSE
#9  157832 b2010010 1285939205      core     1  kalkyl     1 COMPLETED 2010-10-01    TRUE   FALSE
#10 157849 b2010010 1285939205      core     1  kalkyl     1 COMPLETED 2010-10-01    TRUE   FALSE
#..    ...      ...        ...       ...   ...     ...   ...       ...        ...     ...     ...
#Variables not shown: is_multi (lgl), job_type (fctr), duration_sec (int), node_frac (dbl), year
#  (chr), core_hours (dbl), date_block (fctr), proj (chr)

master.prepare.data = function(all.data=FALSE, end.date=study.end.date) {
    logfile = paste0("master.prepare.data_all.data", all.data, "_end.date", end.date, ".log")
    sink(logfile, split = TRUE)
    cat("master.prepare.data: writing log to", logfile, "\n")
    cat("master.prepare.data: all.data =", all.data, "end.date =", end.date, "\n")
    subs = ifelse(all.data, 100, 10)
    do.prepare.data = function(cls, pr, sb, md = study.min.duration) {
        prepare.data(cls, subsample.percent=sb,  min.duration=md, proj=pr, end.date=end.date)  # SNIC bio projects
    }
    d1  = do.prepare.data(c("kalkyl","milou"),           "b",  subs)  # UPPNEX projects
    d2  = do.prepare.data(c("kalkyl","milou","nestor"),  "a",  subs)  # sequencing platform
    d3  = do.prepare.data(c("tintin"),                   "p",  subs)  # p misc projects
    d4  = do.prepare.data(c("tintin"),                   "s",  subs)  # SNIC projects
    d5  = do.prepare.data(c("fysast1"),                  "p",  100)  # p misc projects
    d6  = do.prepare.data(c("halvan"),                   "b",  100)  # UPPNEX projects
    d7  = do.prepare.data(c("halvan"),                   "p",  100)  # p misc projects
    d8  = do.prepare.data(c("kalkyl","milou","tintin","fysast1"), "g", 100)  # g teaching projects
    d9  = do.prepare.data(c("halvan"),                   "bs", 100)  #SNIC Bio projects
    d10 = do.prepare.data(c("kalkyl","tintin"),          "bs", 100)  # SNIC bio projects
    d11 = do.prepare.data(c("fysast1"),                  "s",  100)  # p misc projects
    d12 = do.prepare.data(c("kalkyl"),                   "p",  subs)  # p misc projects
    d13 = do.prepare.data(c("kalkyl"),                   "s",  subs)  # SNIC projects
    sink()
    f = function(x) attr(x, "xdat")
    totaltime.dat <- rbind(f(d1), f(d2), f(d3), f(d4), f(d5), f(d6), f(d7), f(d8), f(d9), f(d10), f(d11), f(d12), f(d13))
    rownames(totaltime.dat) = NULL
    write.data(totaltime.dat, paste0('totaltime_', output.file.tag, '.txt'), gzip = FALSE)
    totaltime.dat <<- totaltime.dat
    bind_rows(d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13) %>% arrange(start_date)
}

count.nodes = function(n, cores, max.cores) {
   # counting nodes in nodes column, number of "," + 1, str_count() from stringr
   # if that's not possible, count nodes via number of cores
   ifelse(n != "", str_count(n, ",") + 1, ceiling(cores / max.cores))
}
make.filenames = function(cluster=default.clusters, sub, dur, proj, data.dir="data", subextension="-COMPLETED", extension=".tsv", gzip=TRUE) {
    if (gzip) extension = paste0(extension, ".gz")
    cluster = match.arg(cluster, valid.clusters)
    if (proj != "") proj = paste0("-", proj) # the way proj is encoded in the filename
    orig.db = paste0(data.dir, "/", cluster, "-orig", proj, subextension, extension)
    prepared.db = paste0(data.dir, "/", cluster, "-processed", proj, subextension, extension)
    prepared.db.filtered = paste0(data.dir, "/", cluster, "-processed", proj, ".sub", sub, ".dur", dur, subextension, extension)
    ####
    invisible(list(orig.db=orig.db, prepared.db=prepared.db, prepared.db.filtered=prepared.db.filtered, gzip=gzip))
}
#for the later files
orig.names = c("date", "job_id", "proj_id", "user", "start", "end", "partition", "cores", "cluster", "nodes", "jobname", "jobstate")
orig.classes = c("NULL","integer","character","NULL","integer","integer","character","numeric","character","character","NULL","character")
#for the 20150806 files
#orig.names = c("job_id","proj_id","user","start","end","partition","cores","cluster","jobstate","nodes")
#orig.classes = c("integer","character","NULL","integer","integer","character","numeric","character","character","character")
orig.rows = 15000000

processed.names = c("job_id","proj_id","start","partition","cores","cluster","jobstate","nodes","start_date","is_core","is_node","is_multi","job_type","duration?sec","node_frac","year","core_hours","date_block")
processed.classes = c("integer","character","integer","character","integer","character","character","integer","character","logical","logical","logical","character","integer","numeric","character","numeric","character")
processed.rows = 6000000

read.orig.data = function(file) {
    dat = tbl_df(read.delim(file, sep="\t", colClasses=orig.classes, header=TRUE, nrows = orig.rows))
    this.cluster = dat$cluster[1]
    max.cores = get.cores(this.cluster)
    n = nrow(dat)
    dat = dat %>% filter(start <= end)
    if (n != nrow(dat)) cat("read.orig.data: removed", n-nrow(dat), "rows for start > end\n")
    dat = dat %>% 
        mutate(start_date = as.Date(as.POSIXct(start, tz="GMT", origin="1970-01-01"))) %>% 
        mutate(nodes = count.nodes(nodes, cores, max.cores)) %>%
        mutate(is_core = (partition %in% c("core","halvan") & cores == 1)) %>%
        mutate(is_node = (nodes == 1 & (partition == "node" | cores == max.cores))) %>%
        mutate(is_multi = (nodes > 1 | cores > max.cores)) %>%
        mutate(job_type = factor(ifelse(is_core, "core", ifelse(is_node, "node", ifelse(is_multi, "multi", "partial"))), levels=valid.jobtype, ordered=TRUE)) %>%
        mutate(duration_sec = end - start) %>%
        mutate(node_frac = round(cores / max.cores, 5)) %>% 
        mutate(year = format(start_date, "%Y")) %>%
        mutate(core_hours = round(cores * duration_sec / 3600., 5)) %>%
        select(-end) %>%
        block.data(blocks.1month) %>%
        arrange(job_id)
}

prepare.data = function(clusters=c("kalkyl","milou"), subsample.percent=10,
                        min.duration=study.min.duration, proj="", end.date=study.end.date,
                        subextension="-COMPLETED") {
    alldat = tbl_df(data.frame())
    xdat = data.frame()
    for (cl in clusters) {
        fn = make.filenames(cl, subsample.percent, min.duration, proj=proj, gzip=TRUE, subextension=subextension)
        cat("prepare.data:", fn$orig.db, "original job data for", cl, proj, "\n")
        dat = read.orig.data(fn$orig.db) %>%
            mutate(proj = proj)
        cat("prepare.data:", fn$orig.db, "nrow =", nrow(dat), "\n")
        if (! is.null(end.date)) {
            end.date = as.Date(as.POSIXct(end.date, tz="GMT"))
            cat("prepare.data:", fn$orig.db, "filtering for end.date <", as.character(end.date), "...\n")
            dat = dat %>%
                filter(as.Date(start_date) < end.date)
            cat("prepare.data:", fn$orig.db, "nrow =", nrow(dat), "\n")
        }
        cat("prepare.data:", fn$orig.db, "writing prepared data to", fn$prepared.db, "...\n")
        write.data(dat, fn$prepared.db, gzip=TRUE)
        sdat = tbl_df(data.frame())
        if (subsample.percent < 100) {
            cat("prepare.data:", fn$orig.db, "with subsampling, excluded jobs are not summarised\n")
            sdat = dat %>% 
                filter(duration_sec >= min.duration) %>% 
                sample_frac(size = subsample.percent/100.) %>% 
                arrange(job_id)
        } else {
            cat("prepare.data:", fn$orig.db, "summarising excluded jobs due to filtering by duration >=",
                min.duration, "\n")
            xxdat = dat %>% filter(duration_sec < min.duration)
            xxdat = xxdat$duration_sec
            xxdat = list(cluster = cl, proj = proj, min.duration = min.duration,
                         n.short = length(xxdat), sec.short = sum(as.numeric(xxdat)))
            cat("prepare.data:", fn$orig.db, "excluded", xxdat$n.short, "too-short jobs totalling",
                xxdat$sec.short, "core sec ==", round(xxdat$sec.short/3600, 3), "core hours, mean",
                round(xxdat$sec.short / xxdat$n.short, 1), "sec / job\n")
            sdat = dat %>% 
                  filter(duration_sec >= min.duration) %>% 
                  arrange(job_id)
            xxdat$n.long = nrow(sdat)
            xxdat$sec.long = sum(as.numeric(sdat$duration_sec))
            xdat = rbind(xdat, xxdat)
        }
        cat("prepare.data:", fn$prepared.db, "duration_sec >=", min.duration, "subsample ",
            subsample.percent,"%, nrow =", nrow(sdat), "\n")
        cat("prepare.data:", fn$prepared.db, "writing filtered data to ", fn$prepared.db.filtered, "...\n")
        write.data(sdat, fn$prepared.db.filtered, gzip=TRUE)
        alldat = rbind(alldat, sdat)
    }
    cat("prepare.data: returning", nrow(alldat), "lines of data for clusters", clusters, "proj",
        ifelse(proj!="", proj, "all"), "\n")
    ####
    alldat = alldat %>% arrange(start_date)
    attr(alldat, "xdat") = xdat
    alldat
}

write.data = function(dat, file, gzip = TRUE) {
    if (gzip) file = gzfile(file, "w")
    write.table(dat, file, sep="\t", col.names=TRUE, row.names=FALSE, quote=TRUE)
    if (gzip) close(file)
}

block.data = function(dat, blocks, column.to.block = "start_date", end.date = study.end.date, start.date) {
    stopifnot(! missing(blocks))
    end.date = as.Date(as.POSIXct(end.date, tz="GMT"))
    if (! missing(start.date)) start.date = as.Date(as.POSIXct(start.date, tz="GMT"))
    blocks = as.Date(as.POSIXct(blocks, tz="GMT"))
    # remove blocks
    new.blocks = blocks[blocks <= end.date]  # for the blocks this must be <= because last must exceed all dates in data
    if (length(new.blocks) < length(blocks)) {
        cat("block.data: removed", length(blocks)-length(new.blocks), "blocks that exceeded", end.date, "\n")
        blocks = new.blocks
    }
    if (! missing(start.date)) {
        new.blocks = blocks[blocks >= start.date]
        if (length(new.blocks) < length(blocks)) {
            cat("block.data: removed", length(blocks)-length(new.blocks), "blocks that exceeded", end.date, "\n")
            blocks = new.blocks
        }
    }
    job_date = as.Date(dat[[column.to.block]])
    # remove jobs but don't be happy about it
    n.jobs_outside = sum(job_date >= end.date)
    if (n.jobs_outside) {
        cat("block.data: WARNING:", n.jobs_outside, "jobs exceed or match", end.date, ", removing...\n")
        dat = dat %>% filter(as.Date(dat[[column.to.block]]) < end.date)
        job_date = as.Date(dat[[column.to.block]])
    }
    if (! missing(start.date)) {
        n.jobs_outside = sum(job_date < start.date)
        if (n.jobs_outside) {
            cat("block.data: WARNING:", n.jobs_outside, "jobs preceded", start.date, ", removing...\n")
            dat = dat %>% filter(as.Date(dat[[column.to.block]]) >= start.date)
            job_date = as.Date(dat[[column.to.block]])
        }
    }
    job_block = character(length(job_date))
    job_block[] = "none"
    for (ib in 1:length(blocks)) {
        in_block = job_date >= blocks[ib] & job_date < blocks[ib + 1]
        job_block[in_block] = as.character(blocks[ib])
    }
    dat$date_block = I(job_block)
    ####
    dat = dat %>%
        filter(date_block != "none") %>%
        mutate(date_block = factor(as.character(date_block), levels=as.character(blocks[-length(blocks)]), ordered=TRUE)) %>% 
        arrange(date_block)
}

