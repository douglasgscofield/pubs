# Functions for loading and integrating probe chl-a fluorescence traces.  Used
# in Richert et al. 
#
# Copyright (c) 2016, Douglas G. Scofield, douglas.scofield@ebc.uu.se
#
# Uppsala University



# load.trace() -- load probe trace files for a given station.
#
# This implementation uses a station.chars dataframe to learn per-station
# naming details.  Loads probe trace files and tweaks the data a bit, in order:
#
# (1) renames columns with inconsistent names
# (2) drops readings below 300 m, depending in 'handle' parameter for station
# (3) drops readings with missing chl-a fluorescence
# (4) drops readings at duplicate depths
# (5) drops readings with outlier chl-a flourescences (>= 1000)

load.trace = function(station) {
    event = station.chars[as.character(station), "event"]
    handle = station.chars[as.character(station), "handle"]
    fn = trace.filename(station, event)
    dat = read.delim(fn)
    # make column names consistent
    cat("load.trace: original names =", names(dat), "\n")
    names(dat)[names(dat) %in% c("Beam.Attenuation", "Beam.Atten..1.m.")] = "BeamAttenuation"
    names(dat)[names(dat) %in% c("FlECO.AFL", "Fluor")] = "chla.fluor"
    names(dat)[names(dat) %in% c("PrDM")] = "depth"
    cat("load.trace: updated names  =", names(dat), "\n")
    if (! "BeamAttenuation" %in% names(dat)) dat$BeamAttenuation = NA
    stopifnot(sum(names(dat) %in% "BeamAttenuation") == 1)
    dat = dat[order(dat$depth), ];
    nr = nrow(dat)
    cat("load.trace:", nr, "rows for station", station, "\n")
    if (handle == "drop300") {
        dat = subset(dat, depth <= 300); nr.new = nrow(dat)
        if (nr > nr.new) cat("load.trace: removed", nr - nr.new,"'drop300' rows for station", station, "\n")
        nr = nr.new
    }
    dat = dat[! is.na(dat$chla.fluor), ]; nr.new = nrow(dat)
    if (nr > nr.new) cat("load.trace: removed", nr - nr.new,"NA rows for station", station, "\n")
    nr = nr.new
    dat = dat[! duplicated(dat$depth), ]; nr.new = nrow(dat)
    if (nr > nr.new) cat("load.trace: removed", nr - nr.new,"duplicated-depth rows for station", station, "\n")
    nr = nr.new
    # check for outlier max: if >10% larger than the next
    dat = subset(dat, abs(dat$chla.fluor) < 1000); nr.new = nrow(dat)
    if (nr > nr.new) cat("load.trace: removed", nr - nr.new,"outlier rows for station", station, "\n")
    nr = nr.new
    # make diffs
    dat$diff = c(0, diff(dat$chla.fluor))
    attr(dat, "station") = station
    attr(dat, "event") = event
    attr(dat, "handle") = handle
    attr(dat, "site.div") = as.list(site.div[station, ])
    #print(head(dat, 3))
    dat
}


# cut.trace() -- integrate probe chl-a fluorescence trace
#
# Uses data loaded by load.trace() to integrate chl-a fluorescence traces.
# Implements various methods of stopping based on stabilisation of the
# integral.
#
# tolerance    : maximum relative difference in integral value allowed over
#                last depth.window m to complete integral
# depth.window : across what depth (m) should be apply tolerance?
#
# See cut.trace2() and its method argument for additional stabilisation
# possibilities.

cut.trace = function(dat, tolerance = 0.01, depth.window = 10) {
    method = "id"
    cat("cut.trace: method =", method, " handle =", attr(dat, "handle"), " tolerance =", tolerance, " depth.window =", depth.window, "\n")
    # adds a 0, 0 line as the first line in dat
    dat = rbind(dat[1, , drop=FALSE], dat, make.row.names = FALSE)
    dat[1,1] = dat[1,2] = 0
    integral = dintegral = numeric(nrow(dat))
    integral[1] = 0
    dintegral[1] = 0
    integral[2] = trapz(dat$depth[1:2], dat$chla.fluor[1:2])
    dintegral[2] = (integral[2] - integral[1]) / (dat[1,"depth"] - dat[2,"depth"])
    for (i in 3:nrow(dat)) {
        integral[i] = trapz(dat$depth[1:i], dat$chla.fluor[1:i])
        dintegral[i] = integral[i] - integral[i - 1]
        # now check to see if integral change > tolerance over last depth.window m
        j = i - 1
        is.stop = 1
        integral.diff = 0
        while (dat$depth[i] < dat$depth[j] + depth.window) {
            delta.m = dat$depth[j+1] - dat$depth[j]
            d.m = (integral[j + 1] - integral[j]) / delta.m
            ########
            msg = sprintf("cut.trace: m= %s  depth= %.1f  i= %d  j,j+1= %d,%d  delta.m= %.1f  %.2f - %.2f d.m= %.4f\n",
                          method, dat$depth[i], i, j, j+1, delta.m, integral[j+1], integral[j], d.m)
            #cat(msg)
            if ((method == "id" && abs(d.m) > tolerance)) {
                is.stop = 0
                break
            }
            j = j - 1
        }
        if (is.stop) {
            cat("cut.trace: is.stop = 1 integral.diff =", integral.diff, " i =", i, " integral[i] =", integral[i], " depth =", dat$depth[i], "\n")
            break
        }
    }
    #cat("cut.trace: making the cut at measurement", i, "at depth", dat$depth[i], "\n")
    mean.at.depth = mean(dat$chla.fluor[(j + 1):i])
    mean.below.depth = mean(dat$chla.fluor[(i + 1):(i + j)])
    cat("cut.trace: mean chla.fluor in depth.window =", mean.at.depth, " below depth.window =", mean.below.depth, "\n")
    dat = dat[1:i, ]
    dat$chla.fluor[2:i] = dat$chla.fluor[2:i] - mean.at.depth
    integral.adjusted = trapz(dat$depth[1:i], dat$chla.fluor[1:i])
    cat("cut.trace: integral =", integral[i], " integral.adjusted =", integral.adjusted, "\n")
    dat$diff = c(0, diff(dat$chla.fluor))
    if (attr(dat, "handle") == "ignore") {
        attr(dat, "integral") = NA
        attr(dat, "integral.adjusted") = NA
    } else {
        attr(dat, "integral") = integral[i]
        attr(dat, "integral.adjusted") = integral.adjusted
    }
    attr(dat, "mean.at.depth") = mean.at.depth
    attr(dat, "mean.below.depth") = mean.below.depth
    #print(head(dat, 3))
    cat("\n")
    dat
}


# cut.trace2() -- integrate probe chl-a fluorescence trace
#
# See cut.trace() for general description.  This method provides additional
# possibilities for method.
#
# method   : stopping method against which tolerance is applied
#   'id'   : relative difference in integrals (method used in cut.trace())
#   'idd'  : relative difference in approximation of integral growth rate
#   'i'    : absolute difference in integrals

cut.trace2 = function(dat, tolerance = 0.01, depth.window = 20, method = c("id", "idd", "i")) {
    method = match.arg(method)
    cat("cut.trace2: method =", method, " handle =", attr(dat, "handle"), " tolerance =", tolerance, " depth.window =", depth.window, "\n")
    # adds a 0, 0 line as the first line in dat
    dat = rbind(dat[1, , drop=FALSE], dat, make.row.names = FALSE)
    dat[1,1] = dat[1,2] = 0
    integral = dintegral = numeric(nrow(dat))
    integral[1] = 0
    dintegral[1] = 0
    integral[2] = trapz(dat$depth[1:2], dat$chla.fluor[1:2])
    dintegral[2] = integral[2] - integral[1]
    for (i in 3:nrow(dat)) {
        # now check to see if integral change > tolerance over last depth.window m
        j = i - 1
        is.stop = 1
        integral.diff = 0
        while (dat$depth[i] < dat$depth[j] + depth.window) {
            ########
            integral.diff = integral[i] - integral[j]
            integral.reldiff = integral.diff / integral[i]
            dintegral[i] = integral.diff
            dintegral.reldiff = dintegral[i] - dintegral[j] / integral.diff
            #dintegral.reldiff = dintegral[i] - dintegral[j] / dintegral[i]
            #dintegral.reldiff = integral.diff
            #dintegral.diff = abs((dintegral[i] - dintegral[j]) / integral[i])
            dintegral.diff = abs((dintegral[i] - dintegral[j]) / (dat$depth[i] - dat$depth[j]))
            msg = sprintf("cut.trace2: m= %s  i,j= %d %d  d= %.1f  %.2f - %.2f = %.4f  i.reldf= %.4f  di.df= %.4f di.reldf= %.4f\n",
                          method, i, j, dat$depth[i], integral[i], integral[j], integral.diff, integral.reldiff, dintegral.diff, dintegral.reldiff)
            cat(msg)
            if ((method == "id" && integral.reldiff > tolerance) ||
                (method == "idd" && dintegral.diff > tolerance) || 
                (method == "i" && integral.diff > tolerance)
                ) {
                #cat("cut.trace2: is.stop = 0 integral.diff =", integral.diff,
                #    " i =", i, " j =", j, " integral[i] =", integral[i],
                #    " depth =", dat$depth[i], "\n")
                is.stop = 0
                break
            }
            j = j - 1
        }
        if (is.stop) {
            cat("cut.trace2: is.stop = 1 integral.diff =", integral.diff, " i =", i, " integral[i] =", integral[i], " depth =", dat$depth[i], "\n")
            break
        }
    }
    #cat("cut.trace2: making the cut at measurement", i, "at depth", dat$depth[i], "\n")
    mean.at.depth = mean(dat$chla.fluor[(j + 1):i])
    mean.below.depth = mean(dat$chla.fluor[(i + 1):(i + j)])
    cat("cut.trace2: mean chla.fluor in depth.window =", mean.at.depth, " below depth.window =", mean.below.depth, "\n")
    dat = dat[1:i, ]
    dat$chla.fluor[2:i] = dat$chla.fluor[2:i] - mean.at.depth
    integral.adjusted = trapz(dat$depth[1:i], dat$chla.fluor[1:i])
    cat("cut.trace2: integral =", integral[i], " integral.adjusted =", integral.adjusted, "\n")
    dat$diff = c(0, diff(dat$chla.fluor))
    if (attr(dat, "handle") == "ignore") {
        attr(dat, "integral") = NA
        attr(dat, "integral.adjusted") = NA
    } else {
        attr(dat, "integral") = integral[i]
        attr(dat, "integral.adjusted") = integral.adjusted
    }
    attr(dat, "mean.at.depth") = mean.at.depth
    attr(dat, "mean.below.depth") = mean.below.depth
    #print(head(dat, 3))
    cat("\n")
    dat
}


# load.data() -- load all probe trace data and integrate with cut.trace()
#
# This implementation requires a focal.stations global variable listing the
# stations.  Returned is a dataframe with integration using depth windows of
# 10, 20 and 30 m (see cut.trace()).

load.data = function(...) {
    sink("load.trace.cut.trace.log", split = TRUE)
    all.data <- list()
    integration.data <- data.frame()
    for (s in focal.stations) {
        cat("\n\nloading data for station", s, "\n\n")
        y <- load.trace(s)
        all.data[[as.character(s)]] <- cut.trace(y, ...)
        x10 <- cut.trace(y, tolerance = 0.01, depth.window = 10)
        x20 <- cut.trace(y, tolerance = 0.01, depth.window = 20)
        x30 <- cut.trace(y, tolerance = 0.01, depth.window = 30)
        r <- list(station = s,
                  event = attr(y, "event"),
                  handle = attr(y, "handle"),
                  ps.10 = attr(x10, "integral"),
                  ps.10.adj = attr(x10, "integral.adjusted"),
                  adj.10 = attr(x10, "mean.at.depth"),
                  ps.20 = attr(x20, "integral"),
                  ps.20.adj = attr(x20, "integral.adjusted"),
                  adj.20 = attr(x20, "mean.at.depth"),
                  ps.30 = attr(x30, "integral"),
                  ps.30.adj = attr(x30, "integral.adjusted"),
                  adj.30 = attr(x30, "mean.at.depth"))
        integration.data <- rbind(integration.data, r)
    }
    rownames(integration.data) = NULL
    attr(all.data, "integration.data") <- integration.data
    sink()
    all.data
}

