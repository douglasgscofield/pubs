# Generate Figure 1 from Dahlö et al, tracking project growth over time.  Uses
# project start-end dates to generate IRanges (from the BioConductor project),
# the 'coverage' of which is used to count active projects.


options(prompt="project_growth> ", stringsAsFactors=FALSE)
#!/sw/apps/R/x86_64/3.0.2/bin/Rscript

#reload = function(doit=FALSE) if (doit) source("project_growth.r")

#source("http://bioconductor.org/biocLite.R")
#biocLite("IRanges")


# rm(list=ls())
library("RSQLite")
library('zoo')
library('IRanges')
library('RColorBrewer')
library('dplyr')
library('stringr')
library('plotrix')
source('legendx.R')

# data file containing names of NGS projects allocated SNIC resources
# does not exist for anonymised data
#bsnic <- scan(file='../data/B_SNIC.txt', what=character())

plot_to_file = 1

cat('
do.plot() to plot figure
')

groups = c('nonNGS', 'NGS','PLATFORM', 'OTHER')

project.type = function(proj_id, g = groups) {
    type = character(length(proj_id))
    type[] = 'OTHER'
    # find all the a projects
    proj.platform = ! is.na(str_match(proj_id, "^a\\d+$")[,1])
    type[proj.platform] = 'PLATFORM'
    # find all the b projects
    proj.uppnex = (! is.na(str_match(proj_id, "^b\\d+$")[,1])) | (proj_id %in% bsnic)
    type[proj.uppnex] = 'NGS'
    # save for general projects
    proj.uppmax = (((! is.na(str_match(proj_id, "^[ps]\\d+$")[,1])) |
                    (! is.na(str_match(proj_id, "^snic.+$")[,1])))    &
                   (! proj_id %in% bsnic))
    type[proj.uppmax] = 'nonNGS'
    stopifnot(all(type %in% groups))
    type
}

do.plot = function() {

    sink('project_growth.log', split=TRUE)

    study.start.year = 2004
    study.start.date = paste0(as.character(study.start.year), "-01-01")
    study.start.numeric = as.numeric(as.POSIXct(study.start.date, tz="GMT"))
    #study.end.date = "2016-07-01"
    study.end.date = "2017-01-01"
    study.end.numeric = as.numeric(as.POSIXct(study.end.date, tz="GMT"))
    window.start.date = "2003-12-01"
    #window.end.date = "2016-07-01"
    window.end.date = "2017-01-01"
    window.start.numeric = as.numeric(as.POSIXct(window.start.date, tz="GMT"))
    window.end.numeric = as.numeric(as.POSIXct(window.end.date, tz="GMT"))
    rate.start.date = "2010-01-01"
    rate.start.numeric = as.numeric(as.POSIXct(rate.start.date, tz="GMT"))

    # connect to the db
    jobs_db = dbConnect(dbDriver("SQLite"), "../data/umax_stats.sqlite")

    # get ngs corehour usage per week
    results <- dbSendQuery(jobs_db,
    "select proj_id,start,end,type from projects;") # create query
    data = fetch(results,n=-1) # get all results
    dbClearResult(results) # reset db

    data = tbl_df(data)
    odata = data

    # convert to epoch time
    data$start_orig = data$start
    data$end_orig = data$end
    data$start = as.numeric(as.POSIXct(strptime(data$start, "%Y%m%d", tz="GMT")))
    data$end = as.numeric(as.POSIXct(strptime(data$end, "%Y%m%d", tz="GMT")))

    axis.start.year = 2006
    axis.end.year = 2017
    total.plot.months = (2017 - axis.start.year) * 12
    cat("total.plot.months", axis.start.year, "to 2016-12 =", total.plot.months, "\n")

    # remove NA rows
    data = data[!is.na(data$start),]

    # select only ngs projects that are newer than 2005 to avoid errors from the projects file
    #data = data[ data$start >= 1041379200 , ]

    #data = subset(data, start >= study.start.numeric)
    data = subset(data, start >= window.start.numeric)

    # select only ngs projects that are started prior to our end date
    #data = subset(data, start < study.end.numeric)
    data = subset(data, start < window.end.numeric)

    # set the end date of all projects ending in the future to end.date
    #data$end[data$end > study.end.numeric] = study.end.numeric
    data$end[data$end > window.end.numeric] = window.end.numeric

    # trim data to only cover the period we have projects
    #first_proj = min(data$start)
    #first_proj = study.start.numeric
    first_proj = window.start.numeric
    data$start = data$start - first_proj
    data$end = data$end - first_proj
    study.end.numeric = study.end.numeric - first_proj

    # convert the data from second resolution to day resolution
    day.sec = 24 * 60 * 60
    data$start = as.integer(data$start / day.sec)
    data$end = as.integer(data$end / day.sec)

    # load list of 'bs' projects, SNIC projects that are NGS-related

    # create a coverage array that is smoothed
    bs_proj = scan("../data/B_SNIC.txt", what = character())
    is.ngs_proj    = data$type == 'uppnex' | data$type == 'platform' | data$proj_id %in% bs_proj
    is.course_proj = data$type == 'course'
    is.staff_proj  = logical(nrow(data))
    staff.proj.pattern = '(test|dummy|staff)'
    is.staff_proj[grep(staff.proj.pattern, data$proj_id)] = TRUE
    is.general_proj    = ! (is.ngs_proj | is.course_proj | is.staff_proj)
    cat("total general projects =", sum(is.general_proj), "\n")
    cat("total ngs projects =", sum(is.ngs_proj), "\n")
    cat("total staff projects", staff.proj.pattern, "=", sum(is.staff_proj), "\n")
    cat("total course projects =", sum(is.course_proj), "\n")
    c.general = unlist(as.list(coverage(IRanges(start = data$start[is.general_proj],
                                                end   = data$end[is.general_proj]))))
    c.ngs     = unlist(as.list(coverage(IRanges(start = data$start[is.ngs_proj],
                                                end   = data$end[is.ngs_proj]))))
    c.course  = unlist(as.list(coverage(IRanges(start = data$start[is.course_proj],
                                                end   = data$end[is.course_proj]))))

    crossover.day.n = sum((c.ngs - c.general) < 0) + 1
    cat("general < ngs crossover.day =", crossover.day.n, " "); print(as.Date(window.start.date) + crossover.day.n)
    cat("crossover.n general:", c.general[crossover.day.n], " ngs:", c.ngs[crossover.day.n], "\n")
    cat("max general n =", max(c.general), "\n")
    cat("max ngs n =", max(c.ngs), "\n")
    cat("max course n =", max(c.course), "\n")

    first.ngs.project.day = min(which(c.ngs >= 1))
    cat("first.ngs.project.day =", first.ngs.project.day, " "); print(as.Date(window.start.date) + first.ngs.project.day)
    ten.ngs.project.day = min(which(c.ngs >= 10))
    cat("ten.ngs.project.day =", ten.ngs.project.day, " "); print(as.Date(window.start.date) + ten.ngs.project.day)
    fifty.ngs.project.day = min(which(c.ngs >= 50))
    cat("fifty.ngs.project.day =", fifty.ngs.project.day, " "); print(as.Date(window.start.date) + fifty.ngs.project.day)
    fifty.ngs.project.day = min(which(c.ngs >= 50))
    cat("fifty.ngs.project.day =", fifty.ngs.project.day, " "); print(as.Date(window.start.date) + fifty.ngs.project.day)

    k = 30  # right-aligned rolling mean window size

    general_smooth = rollmean(c.general, k, align="right")
    ngs_smooth = rollmean(c.ngs, k, align="right")
    course_smooth = rollmean(c.course, k, align="right")
    plotdat = data.frame(x = 1:length(general_smooth), general_smooth, ngs_smooth, course_smooth)
    days.in.year = function(y) {
        dim = function(y, m) {  # calculate the number of days in month m of year y
            d = c(31, NA, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
            if (! is.na(d[m])) return(d[m])
            if      (y %% 400 == 0) return(29)
            else if (y %% 100 == 0) return(28)
            else if (y %% 4 == 0)   return(29)
            else                    return(28)
        }
        diy = function(y) sum(sapply(1:12, function(m) dim(y, m)))
        diy(y)
    }

    # compute columns of tick positions and tick labels for plotting
    plotdat$ticks = plotdat$subticks = NA
    plotdat$labels = ""
    iy = 1
    start.iy = 1
    for (y in 2004:2017) {
        plotdat$ticks[iy] = iy
        plotdat$labels[iy] = as.character(y)
        if (y == axis.start.year)
            start.iy = iy
        if (y < axis.end.year) {
            nd = days.in.year(y)
            subnd = round((1:3) * (nd / 4))
            plotdat$subticks[iy + subnd] = iy + subnd
            iy = iy + nd
        }
    }
    plotdat = plotdat[start.iy:nrow(plotdat), ]

    rawdat = data.frame(c.general, c.ngs, c.course)
    rawdat = rawdat[(start.iy - k + 1):nrow(rawdat), ]

    stats = list()
    for (p in c('general', 'ngs', 'course')) {
        rawnm = paste0('c.', p)
        snm = paste0(p, '.num_at_start'); enm = paste0(p, '.num_at_end'); rnm = paste0(p, '.monthly_rate')
        stats[[snm]] = rawdat[1, rawnm]
        stats[[enm]] = rawdat[nrow(rawdat), rawnm]
        stats[[rnm]] = (stats[[enm]] - stats[[snm]]) / total.plot.months
        dnm.1 = paste0(p, '.n1.day')
        snm.1 = paste0(p, '.n1.num_at_start')
        rnm.1 = paste0(p, '.n1_monthly_rate')
        stats[[dnm.1]] = min(which(rawdat[[rawnm]] >= 1))
        stats[[snm.1]] = rawdat[stats[[dnm.1]], rawnm]
        stats[[rnm.1]] = (stats[[enm]] - stats[[snm.1]]) / (total.plot.months - ((stats[[dnm.1]] - 1) * 12 / 365.25))
        dnm.10 = paste0(p, '.n10.day')
        snm.10 = paste0(p, '.n10.num_at_start')
        rnm.10 = paste0(p, '.n10_monthly_rate')
        stats[[dnm.10]] = min(which(rawdat[[rawnm]] >= 10))
        stats[[snm.10]] = rawdat[stats[[dnm.10]], rawnm]
        stats[[rnm.10]] = (stats[[enm]] - stats[[snm.10]]) / (total.plot.months - ((stats[[dnm.10]] - 1) * 12 / 365.25))
        if (p != 'course') {
            dnm.50 = paste0(p, '.n50.day')
            snm.50 = paste0(p, '.n50.num_at_start')
            rnm.50 = paste0(p, '.n50_monthly_rate')
            stats[[dnm.50]] = min(which(rawdat[[rawnm]] >= 50))
            stats[[snm.50]] = rawdat[stats[[dnm.50]], rawnm]
            stats[[rnm.50]] = (stats[[enm]] - stats[[snm.50]]) / (total.plot.months - ((stats[[dnm.50]] - 1) * 12 / 365.25))
        }
        smoothnm = paste0(p, '_smooth')
        snm = paste0('smooth.', snm); enm = paste0('smooth.', enm); rnm = paste0('smooth.', rnm)
        stats[[snm]] = plotdat[1, smoothnm]
        stats[[enm]] = plotdat[nrow(plotdat), smoothnm]
        stats[[rnm]] = (stats[[enm]] - stats[[snm]]) / total.plot.months
    }
    sapply(names(stats), function(.x) cat(.x, "=", stats[[.x]], "\n"))

    np = as.Date(rate.start.date) - as.Date(window.start.date) + 1
    cat("NGS projects at", rate.start.date, "=", c.ngs[np], "\n")
    cat("General projects at", rate.start.date, "=", c.general[np], "\n")
    cat("Course projects at", rate.start.date, "=", c.course[np], "\n")
    j = (2016 - 2010 + 1) * 12 # months from 2010-01-01 to 2017-01-01
    cat("NGS project growth rate from", rate.start.date, "=", (stats$ngs.num_at_end - c.ngs[np]) / j, "\n")
    cat("General project growth rate from", rate.start.date, "=", (stats$general.num_at_end - c.general[np]) / j, "\n")
    cat("Course project growth rate from", rate.start.date, "=", (stats$course.num_at_end - c.course[np]) / j, "\n")


    cat("first plot points\n")
    print(head(plotdat, 5))
    cat("last plot points\n")
    print(tail(plotdat, 5))

    general_color = brewer.pal(9, "Reds")[7] #'#123763'
    ngs_color = brewer.pal(9, "Greens")[7] #'#cd0000'
    general_color_trans = paste0(general_color, "33")
    ngs_color_trans = paste0(ngs_color, "33")
    general_lty = 2
    ngs_lty = 1


    filename = 'Fig_projects_stacked.pdf'
    resx = 85 / 25.4
    resy = 90 / 25.4
    factor = resx / 10
    if(plot_to_file == 1){
        pdf(filename, width=resx, height=resy, pointsize=10)
    }

    # initiate plot
    if (0) {
        opa = par(mar=c(2.0, 2.7, 0.5, 0.4), lwd=2, mgp=c(1.8, 0.4, 0),
                  xaxs='i', yaxs='i', las=1, ps=9, tcl=-0.2, cex.axis=0.8)
    } else if (0) {
        ## opa = par(mfcol=c(3, 1), mar=c(2.0, 2.7, 0.5, 0.4), lwd=2, mgp=c(1.8, 0.4, 0),
        ##           xaxs='i', yaxs='i', las=1, ps=9, tcl=-0.2, cex.axis=0.8)
        layout(matrix(c(1, 2, 3), byrow=TRUE, nrow=1, ncol=3),
               widths=c(12, 8, 5))
        opa = par(mar=c(2.0, 2.7, 0.5, 0.4), lwd=2, mgp=c(1.8, 0.4, 0),
                  xaxs='i', yaxs='i', las=1, ps=9, tcl=-0.2, cex.axis=0.8)
    } else {
        layout(matrix(c(1, 1, 1, 0, 2, 2, 0, 0, 3), byrow=TRUE, nrow=3, ncol=3),
               widths=c(128, 96, 168))
        opa = par(mar=c(2.0, 2.5, 0.5, 0.4), oma=c(0, 1, 0, 0.5),
                  lwd=2, mgp=c(1.8, 0.4, 0),
                  xaxs='i', yaxs='i', las=1, ps=9, tcl=-0.2, cex.axis=0.8)
    }

    #
    # project numbers
    #

    this.xlim = range(plotdat$x)
    step = 100
    this.ylim = c(0, ceiling(max(plotdat$general_smooth, plotdat$ngs_smooth) / step) * step)
    plot(1, type='n', xlim=this.xlim, ylim=this.ylim+c(0,10), axes=F, xlab="", ylab="")

    n = nrow(plotdat)

    # calculate Y axis
    # Y axis
    yat = seq(this.ylim[1], this.ylim[2], by=step)

    # calculate X axis
    # X axis
    ### generate the year labels
    # get the day of the year for the first project
    doj = as.integer(strftime(as.POSIXct(first_proj, origin="1970-01-01"), format = "%j"))
    # get the first year number to plot
    first_year = as.integer(strftime(as.POSIXct(first_proj, origin="1970-01-01"), format = "%Y")) + 1

    # create an X-axis with 1 year longer in each direction than fits in the plot
    xat = seq(-365.25,length(plotdat$general_smooth)+365.25, 365.25) + (365.25 - doj)

    # grid lines
    abline( h = yat, lty = 3, lwd=1.0, col = brewer.pal(9, "Greys")[4] )

    # plot the ngs area under curve
    uplotdat = subset(plotdat, ngs_smooth > 0)
    un = nrow(uplotdat)

    # plot the gen line
    lines(plotdat$x, plotdat$general_smooth, col=general_color, lty=general_lty)

    # plot the ngs line
    lines(uplotdat$x, uplotdat$ngs_smooth, col=ngs_color, lty=ngs_lty)

    # axes
    labels = c(NA, seq(first_year, first_year+length(xat)-2))
    labels[length(labels) - 1] = NA  # blank 2017
    labels[length(labels)] = NA  # blank 2018
    xi = 2*(1:length(labels)) - 1  # label positions in new axis
    xlabels = rep(NA, 2*length(labels) - 1)
    xlabels[xi] = labels

    xat.which = which(! is.na(plotdat$ticks))
    xat.ticks = plotdat$ticks[xat.which]
    xat.labels = plotdat$labels[xat.which]
    xat.subticks = plotdat$subticks[! is.na(plotdat$subticks)]
    axis(side=1, at=xat.subticks, labels=FALSE, tcl=-0.1, lwd.ticks=0.5, lwd=0)
    axis(side=1, at=xat.ticks, labels=FALSE)
    text(xat.ticks, 1, labels=xat.labels, srt=45, adj=c(1.3, 1.5), xpd=NA, cex=0.8)

    ylabels=yat
    axis(side=2, at=yat, labels=c(0, ylabels[2:length(ylabels)]), cex=0.7)

    # header
    title(ylab="Active projects", cex.lab=1.2, xpd=NA)

    # legend
    legendx('topleft', c('NGS projects', 'Non-NGS projects'),
           lwd=2, col=c(ngs_color, general_color), seg.len=3, lty=c(ngs_lty, general_lty),
           bty="o", pch=NA, cex=1.0, box.col="white", bg="white", inset=c(0.10, -0.15), box.hex=0.8, xpd=NA)
    usr = par('usr')
    boxed.labels(usr[1]-0.1*(usr[2]-usr[1]),usr[4]-0.02*(usr[4]-usr[3]), "A",
                 cex=1.1, font=2, bg="white", border=NA, xpd=NA)

    #
    # unique PIs
    #
    opa.b=par(mar=c(2.0, 2.5, 1.0, 0.4))

    uppm = unique_pis_per_month
    m = substring(as.character(uppm$month), 1, 4)
    this.ylim = range(c(0, round(uppm[,2:3]+70, 2)))
    plot(1, 1, type='n', xlim=c(1, length(m)), ylim=this.ylim, axes=F, xlab="", ylab="")
    #plot(1, 1, type='n', xlim=c(1, length(m)), ylim=range(c(0, uppm[,2:3])), axes=F, xlab="", ylab="")
    # axes
    xat.ticks = seq(1, length(m), by=12)
    xat.subticks = seq(1, length(m), by=3); xat.subticks = xat.subticks[! xat.subticks %in% xat.ticks]
    axis(1, at=xat.ticks,    labels=FALSE, lwd=1)
    axis(1, at=xat.subticks, labels=FALSE, tcl=-0.1, lwd.ticks=0.5, lwd=0)
    text(xat.ticks, 1, labels=m[xat.ticks], srt=45, adj=c(1.3, 1.5), xpd=NA, cex=0.8)
    yat = axTicks(2)
    #title(ylab=expression("Total storage (Petabytes "==10^{15}*" bytes)"), line=1.5)
    title(ylab="Principal investigators", line=1.8, cex.lab=1.2, xpd=NA)
    # grid lines
    abline( h = yat[-1], lty = 3, lwd=1.0, col = brewer.pal(9, "Greys")[4] )
    # y-axis
    axis(2, at = yat)
    # storage
    lines(seq_along(m), uppm$uniqpi.NGS, lwd=2, col=ngs_color, lty=ngs_lty)
    lines(seq_along(m), uppm$uniqpi.nonNGS, lwd=2, col=general_color, lty=general_lty)
    usr = par('usr')
    boxed.labels(usr[1]-0.72*(usr[2]-usr[1]),usr[4]-0.02*(usr[4]-usr[3]), "B",
                 cex=1.1, font=2, bg="white", border=NA, xpd=NA)
    par(opa.b)

    #
    # total storage
    #
    opa.c=par(mar=c(2.0, 2.5, 1.0, 0.4))

    m = sort(c(substring(as.character(monthly$date), 1, 7), '2013-11'))
    xm = (1:length(m))[m != '2013-11']
    m = substring(m, 1, 4)
    this.ylim = range(0, round((max(monthly[,2:3])+256)/512, 0)*512)
    plot(1, 1, type='n', xlim=c(1, length(m)+1), ylim=this.ylim+c(0, 50), axes=F, xlab="", ylab="")
    #plot(1, 1, type='n', xlim=c(1, length(m)+1), ylim=range(monthly[,2:5]), log='y', axes=F, xlab="", ylab="")
    # axes
    xat.ticks = seq(1, length(m)+1, by=12)
    xat.subticks = seq(1, length(m), by=3); xat.subticks = xat.subticks[! xat.subticks %in% xat.ticks]
    axis(1, at=xat.ticks,    labels=FALSE, lwd=1)
    axis(1, at=xat.subticks, labels=FALSE, tcl=-0.1, lwd.ticks=0.5, lwd=0)
    text(xat.ticks, 1, labels=c(m,'2017')[xat.ticks], srt=45, adj=c(1.3, 1.5), xpd=NA, cex=0.8)
    #yat = axis(2, labels=FALSE)
    #axis(2, at=yat, labels=sprintf("%.1f", yat/1024), lwd=NA)
    yat = seq(this.ylim[1], this.ylim[2], by=512)
    cat('yat =', yat, '\n')
    #title(ylab=expression("Total storage (Petabytes "==10^{15}*" bytes)"), line=1.5)
    title(ylab="Total storage (PiB)", line=1.7, cex.lab=1.2, xpd=NA)
    # grid lines
    abline( h = yat[-1], lty = 3, lwd=1.0, col = brewer.pal(9, "Greys")[4] )
    # y-axis
    axis(2, at=yat, labels=sprintf("%.1f", yat/1024))
    # storage
    lines(xm, monthly$NGS, lwd=2, col=ngs_color, lty=ngs_lty)
    lines(xm, monthly$nonNGS, lwd=2, col=general_color, lty=general_lty)
    # tag
    usr = par('usr')
    boxed.labels(usr[1]-2.0*(usr[2]-usr[1]),usr[4]-0.02*(usr[4]-usr[3]), "C",
                 cex=1.1, font=2, bg="white", border=NA, xpd=NA)
    par(opa.c)


    if(plot_to_file == 1){
        # write to file
        dev.off()
    }

    par(opa)

    sink()
}



do.storage <- function() {

    dateStart = "2010-01-01"
    dateEnd = "2017-01-01"

    # using the ram drive
    #general_db = dbConnect(dbDriver("SQLite"),"/dev/shm/dahlo/general.sqlite")
    general_db = dbConnect(dbDriver("SQLite"),"../data/umax_stats.sqlite")

    # the storage usage for all projects
    results <- dbSendQuery(general_db, paste("SELECT * FROM storage_history WHERE date>='", dateStart, "' AND date<'", dateEnd, "'", sep='')) # create query

    storage_table = fetch(results,n=-1) # get all results
    dbClearResult(results) # reset db
    storage_table = tbl_df(storage_table)
    cat('do.storage :', nrow(storage_table), 'projects between', dateStart, 'and', dateEnd, '\n')

    storage_table$type = project.type(storage_table$proj_id)

    # divide the data into groups
    storage_usage = list()
    storage_groups = split(storage_table, storage_table$type)

    stopifnot(all(names(storage_groups) %in% groups) && all(groups %in% names(storage_groups)))

    # for each group
    for(group in names(storage_groups)){
        current_group = with(storage_groups[[group]], tapply(usage+noback_usage, date, FUN=sum, na.rm=TRUE))
        current_group = data.frame(as.Date(names(current_group)), current_group)
        colnames(current_group) = c('date', group)
        rownames(current_group) = c()
        storage_usage[[group]] = current_group
    }

    # merge all the groups into a single data frame
    master = data.frame(date=as.Date(character()), usage=numeric())
    for(group in names(storage_groups)){
        master = merge(master, storage_usage[[group]], by='date', all=TRUE)
    }
    # remove the first column, it's filled with only NAs
    master = master[,-2]

    #  construct month identifier
    master$month = format(master$date,format="%y %B")

    # compute monthly averages for all groups
    monthly = data.frame(date=rep(NA, length(unique(master$month))))
    for(group in groups){
        current_average = tapply(master[,group], master$month, FUN=mean, na.rm=TRUE)
        monthly[,'date'] = names(current_average)
        monthly[, group] = current_average
    }

    # convert back to dates
    monthly$date = as.Date(paste('01', monthly$date, sep=' '), format="%d %y %B")
    monthly = monthly[order(monthly$date),]

    # save.image(file='saved_ws_after_script_run_170222.RData')
    cat('do.storage :', nrow(monthly), 'monthly results returned\n')
    invisible(monthly)
}

do.pis <- function()
{
    general_db = dbConnect(dbDriver("SQLite"), "../data/umax_stats.sqlite")
    results <- dbSendQuery(general_db, paste("SELECT proj_id,pi,start,end FROM projects")) # create query
    pis = fetch(results,n=-1) # get all results
    dbClearResult(results) # reset db
    pis = tbl_df(pis)

    pis = pis %>%
          mutate(type = project.type(proj_id)) %>%
          filter(type != 'OTHER') %>%
          filter(! is.na(pi), pi != "") %>%
          mutate(start = as.Date(start, format="%Y%m%d"),
                 end = as.Date(end, format="%Y%m%d")) %>%
          filter(! is.na(start), ! is.na(end))

    pi.uni = unique(pis$pi)
    pi.map = setNames(as.integer(str_match(pi.uni, "\\((\\d+)\\)")[,2]),
                      pi.uni)
    pi.id = !is.na(pi.map)
    cat('do.pis : pi.map first scan,', length(pi.map), 'PIs,', sum(!pi.id),'have no IDs\n')

    # fuzzy match
    pi.dist = adist(pi.uni[pi.id],pi.uni[!pi.id], partial = TRUE, ignore.case = TRUE)
    matches = which(pi.dist<=7, arr.ind=TRUE)
    pi.map[pi.uni[!pi.id][matches[, 2]]] = pi.map[pi.uni[pi.id][matches[, 1]]]
    cat('do.pis : pi.map after fuzzy match,', sum(is.na(pi.map)),'have no IDs\n')

    # for the ones that still are not found, a manual curation was done
    curated = read.table(file='curated.txt', col.names=c("pi", "id"), sep="\t", stringsAsFactors=FALSE)
    curated.id = !is.na(curated[,2])
    curated[!curated.id, 2] = sample(seq(11000, length.out=sum(!curated.id)))
    pi.map[curated[, 1]] = curated[, 2]
    cat('do.pis : pi.map after applying curated.txt,', sum(is.na(pi.map)),'have no IDs\n')

    m = seq(as.Date("2010-01-01"), as.Date("2017-01-01"), by="month")
    ans = data.frame(month=m,
                     uniqpi.nonNGS=0,      uniqpi.NGS=0,       # number of unique PIs
                     projpi.nonNGS=0,      projpi.NGS=0,       # mean projects/PI
                     oneprojpi.nonNGS=0,   oneprojpi.NGS=0,    # number of PIs with 1 project
                     multiprojpi.nonNGS=0, multiprojpi.NGS=0)  # mean projects/PI for PIs with >1 project
    for(i in 1:length(m)){
        # get all projects overlapping the current month
        current_month.nonNGS = subset(pis, start <= m[i] & end >= m[i] & type == 'nonNGS')
        current_month.NGS = subset(pis, start <= m[i] & end >= m[i] & type == 'NGS')
        projcount.nonNGS = table(pi.map[current_month.nonNGS$pi])
        projcount.NGS = table(pi.map[current_month.NGS$pi])
        ans$projpi.nonNGS[i] = mean(projcount.nonNGS)
        ans$projpi.NGS[i] = mean(projcount.NGS)
        ans$uniqpi.nonNGS[i] = length(unique(pi.map[current_month.nonNGS$pi]))
        ans$uniqpi.NGS[i] = length(unique(pi.map[current_month.NGS$pi]))
        ans$oneprojpi.nonNGS[i] = sum(projcount.nonNGS == 1)
        ans$oneprojpi.NGS[i] = sum(projcount.NGS == 1)
        ans$multiprojpi.nonNGS[i] = mean(projcount.nonNGS[projcount.nonNGS > 1])
        ans$multiprojpi.NGS[i] = mean(projcount.NGS[projcount.NGS > 1])
        if (i == length(m)) {
            cat('mean +- SE for other PI > 1 project =',  mean(projcount.nonNGS[projcount.nonNGS > 1]), "+-",
                sqrt(var(projcount.nonNGS[projcount.nonNGS > 1])/sum(projcount.nonNGS > 1)), "\n")
            cat('mean +- SE for NGS PI > 1 project =',  mean(projcount.NGS[projcount.NGS > 1]), "+-",
                sqrt(var(projcount.NGS[projcount.NGS > 1])/sum(projcount.NGS > 1)), "\n")
            print(wilcox.test(projcount.nonNGS[projcount.nonNGS > 1], projcount.NGS[projcount.NGS > 1], conf.int=TRUE, exact=FALSE))
        }
    }
    invisible(ans)
}

if (readline('Should the monthly dataframe be recalculated? [Y/n] ') == 'Y') {
    monthly <- do.storage()
} else cat('Not recalculating')

if (readline('Should the unique_pis_per_month dataframe be recalculated? [Y/n] ') == 'Y') {
    unique_pis_per_month <- do.pis()
} else cat('Not recalculating')

