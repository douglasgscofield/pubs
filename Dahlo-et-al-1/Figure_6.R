# $ module load R/3.2.3

# example run:
# Rscript Figure_6.R '2016-01-01' '2017-01-01'


# rm(list=ls())
# library('RMySQL')
library("RSQLite")
library('stringr')
library('foreach')
library('doMC')
library('RColorBrewer')
library('MASS')
library('hash')
library('KernSmooth')
library('abkde')




# for devel use
dateStart = "2016-01-01"
dateEnd = "2017-01-01"

# get arguments
args = args<-commandArgs(TRUE)
dateStart = args[1]
dateEnd = args[2]


# connect to the db
jobs_db = dbConnect(dbDriver("SQLite"),"../data/efficiency.sqlite")




# # get recent milou jobs
results <- dbSendQuery(jobs_db, paste("select * from jobs where cluster='milou' and date_finished>='", dateStart, "' and date_finished<'", dateEnd, "'", sep='')) # create query

jobs = fetch(results,n=-1) # get all results
dbClearResult(results) # reset db


# add jobs from tintin
results <- dbSendQuery(jobs_db, paste("select * from jobs where cluster='tintin' and date_finished>='", dateStart, "' and date_finished<'", dateEnd, "'", sep='')) # create query
jobsAdd = fetch(results,n=-1) # get all results
dbClearResult(results) # reset db
jobs = rbind(jobs,jobsAdd)


# add jobs from nestor
results <- dbSendQuery(jobs_db, paste("select * from jobs where cluster='nestor' and date_finished>='", dateStart, "' and date_finished<'", dateEnd, "'", sep='')) # create query
jobsAdd = fetch(results,n=-1) # get all results
dbClearResult(results) # reset db
jobs = rbind(jobs,jobsAdd)


# add jobs from halvan
results <- dbSendQuery(jobs_db, paste("select * from jobs where cluster='halvan' and date_finished>='", dateStart, "' and date_finished<'", dateEnd, "'", sep='')) # create query
jobsAdd = fetch(results,n=-1) # get all results
dbClearResult(results) # reset db
jobs = rbind(jobs,jobsAdd)




# skip jobs with only 1 count
jobs = jobs[jobs$counts>1,]



### condense the data frame to one proj per row

# calculate the size of the jobs
jobs$corecounts = jobs$cores * jobs$counts

# summarize everything per project
mem = with(jobs, tapply(mem_peak/mem_limit *100 *corecounts, proj_id, FUN = sum))
data = data.frame(pid = names(mem), mem = mem, cpu = with(jobs, tapply(cpu_mean *corecounts, proj_id, FUN = sum)), corecounts = with(jobs, tapply(corecounts, proj_id, FUN = sum)), counts = with(jobs, tapply(counts, proj_id, FUN = sum)), stringsAsFactors = F, row.names=NULL)


# calculate the average again
data$mem = data$mem/data$corecounts
data$cpu = data$cpu/data$corecounts



# define the group names and order
groups = c('UPPMAX', 'UPPNEX','PLATFORM')

# find all the a projects
lv = is.na(str_match(data$pid, "^plat\\d+$"))
data$type[!lv] = groups[3]

# find all the b projects
lv = is.na(str_match(data$pid, "^bio\\d+$"))
data$type[!lv] = groups[2]

# find all the uppmax projects, ie not a or b projs
lv = is.na(str_match(data$pid, "^gen\\d+$"))
data$type[!lv] = groups[1]



# used to type all the jobs, for manual data mucking
if(0){

    # find all the a projects
    lv = is.na(str_match(jobs$proj_id, "^plat\\d+$"))
    jobs$type[!lv] = groups[3]

    # find all the b projects
    lv = is.na(str_match(jobs$proj_id, "^bio\\d+$"))
    jobs$type[!lv] = groups[2]

    # find all the uppmax projects, ie not a or b projs
    lv = is.na(str_match(jobs$proj_id, "^gen\\d+$"))
    jobs$type[!lv] = groups[1]


    # prepare data
    jobs$core_usage = jobs$cpu_mean * jobs$cores
    jobs$cpu_efficiency = jobs$core_usage / (jobs$cores*100)
    jobs$mem_efficiency = jobs$mem_peak / jobs$mem_limit
    jobs$job_efficiency = pmax(jobs$cpu_efficiency, jobs$mem_efficiency)



    ### numbers used in the paper under Results and Discussion - Efficiency

    # cpu mean/median
    # mean(data$cpu[data$type=='UPPNEX'])
    # mean(data$cpu[data$type=='UPPMAX'])
    quantile(data$cpu[data$type=='UPPNEX'])
    quantile(data$cpu[data$type=='UPPMAX'])
    quantile(data$cpu[data$type=='PLATFORM'])
    # quantile(jobs$cpu_efficiency[jobs$type=='UPPNEX'])
    # quantile(jobs$cpu_efficiency[jobs$type=='UPPMAX'])
    # quantile(jobs$cpu_efficiency[jobs$type=='PLATFORM'])

    # memory mean/median
    # mean(data$mem[data$type=='UPPNEX'])
    # mean(data$mem[data$type=='UPPMAX'])
    quantile(data$mem[data$type=='UPPNEX'])
    quantile(data$mem[data$type=='UPPMAX'])
    quantile(data$mem[data$type=='PLATFORM'])

    




    # mean(jobs$mem_peak[jobs$type=='UPPNEX' & jobs$state=='COMPLETED'])
    # median(jobs$mem_peak[jobs$type=='UPPNEX' & jobs$state=='COMPLETED'])
    # mean(jobs$mem_peak[jobs$type=='UPPMAX' & jobs$state=='COMPLETED'])
    # median(jobs$mem_peak[jobs$type=='UPPMAX' & jobs$state=='COMPLETED'])

}


cpu = weighted_data$cpu
mem = weighted_data$mem



### create the surface plot

# bins
nbins <- 2000 # production

# KernSmooth kde with abkde estimated bandwidth
surf <- lapply(split(data.frame(cpu, mem), weighted_data$type),
    function(x){

            density <- fbkde(x, .verbose=FALSE) # Estimate bandwidth by maximum likelihood
            bandwidth <- density$a[which.max(density$w)] # save optimal bandwidth
            bkde2D(x, bandwidth=bandwidth, gridsize=c(nbins,nbins), list(c(0,100),c(0,100)), truncate = TRUE )$fhat
        })


ax = seq(0 , 100, length.out=nbins+1)
ay = seq(0 , 100, length.out=nbins+1)

# for both uppmax and uppnex
# get order from names(split(data.frame(cpu, mem), weighted_data$type))
names = c('Non-NGS projects', 'NGS projects', 'Platform projects')
# colors = c("Greens", 'Red')
colors = c('YlOrBr', 'BuGn', 'Blues')
plot_order = c("UPPNEX", "UPPMAX")






# prepare to scale the colorRampPalette after the size of category.
count_sum = rep(NA, length(surf))
for(i in 1:length(surf)){
	count_sum[i] = sum(data$corecount[data$type==groups[i]])

}














# color dots by analysis type?
classify = 0
weighted_data$class_col = "black"
if(classify == 1){

    # read the file
    class = readLines("../data/B_SNIC.txt")

    # init
    types = c(1,4:8) # De novo, RNA, Methylation, WGS, Metagenomics
    type_col = brewer.pal(length(types),"Accent")

    # parse each line
    for(line in class){

        # remove all spaces
        line = gsub(" ", "", line)

        # split on tabs
        line = unlist(strsplit(line, "\t"))


        # process each analysis type
        for(type in 1:length(types)){

            # split the project list on commas
            current = unlist(strsplit(line[types[type]], ','))

            # for each project
            for(proj in current){

                # check if the project is present in the statistics
                weighted_data$class_col[weighted_data$pid==proj] = type_col[type]

            }
        }
    }
}
















# customize plot
# color gradient
plot_opt_ncol = 10 # the number of colors the heat map should have
plot_opt_image = 1 # plot the heat map?
plot_opt_scale_image = 0 # scale the heat man intensity after group size?
plot_opt_points = 1 # plot the data points?
plot_opt_subplot_titles = 1 # plot title for each subplot?
plot_opt_contour = 0 # plot the contours in the heat map?
plot_opt_grid = 1 # plot a grid in the plot?
plot_opt_axis = 1 # plot axes in the plot?
plot_opt_save_to_file = 1 # save to a file?# plot uppmax and uppnex
plot_opt_file_name = 'Figure_6_efficiency.pdf'



if(plot_opt_save_to_file == 1){
    # save to file
    pdf(plot_opt_file_name, 4.33, 2.55)
    cex = 0.7
    # set the plot properties, 2 sub-plots, margins around in outer and inner area
    par(mfrow=c(1,length(plot_order)), las=1, mar=c(2, 2, 2, .5), oma=c(1, 1.25, 0, .25), cex=cex)
}else{
    cex=2.2
    # set the plot properties, 2 sub-plots, margins around in outer and inner area
    par(mfrow=c(1,length(plot_order)), las=1, mar=c(2, 3, 3, .5), oma=c(3, 2, 1, 1), cex=cex)
}

for(i in 1:length(plot_order)){

	# pick out the index of surf that contains the project type we are interested in
    surf_idx = which(names(surf)==plot_order[i])
    g_idx = which(groups==plot_order[i])


    if(plot_opt_image == 1){

    	# generate colors
    	col = colorRampPalette(c("white", brewer.pal(9, colors[g_idx])))(plot_opt_ncol)

    	if(plot_opt_scale_image == 1){
	    	# scale the colorRampPalette after the size of category.
	    	col = col[1: ceiling(plot_opt_ncol * count_sum[surf_idx] / max(count_sum)) ]
	    }


        # plot the surface
        image(ax, ay, surf[[surf_idx]], col=col , useRaster=TRUE, xlab="", ylab="", axes=FALSE)

    }

    if(plot_opt_points == 1){
        # plot data points
        points(cpu[weighted_data$type == groups[g_idx]], mem[weighted_data$type == groups[g_idx]], pch=16, cex=.35, col=weighted_data$class_col[weighted_data$type == groups[g_idx]])
    }

    if(plot_opt_subplot_titles == 1){
        # add title to subplot
        title(names[g_idx])
    }


    if(plot_opt_contour == 1){
        # draw a contour
        contour(x=ax[1:(length(ax)-1)], y=ay[1:(length(ay)-1)], z=surf[[surf_idx]], add=T, lwd=.3, drawlabels=FALSE, nlevels=10)
    }


    # create axis marks
    markers_x = seq(0,100,25)
    markers_y = seq(0,100,25)

    if(plot_opt_grid == 1){
        # print grid lines
        abline( v = seq(0,100,25), lty = 3, col = brewer.pal(9, "Greys")[4], lwd=.3 )
        abline( h = seq(0,100,25), lty = 3, col = brewer.pal(9, "Greys")[4], lwd=.3 )
    }


    if(plot_opt_axis == 1){
        # draw logarithmic x axis (subtract 0.005 to get the last mark in the plot)
        axis(1, at=markers_x, lab=markers_x, tck=-.025, cex.axis=0.95, mgp = c(3, 0.5, 0))
        # axis(1, at=c(0, markers_x, lab=c(0,markers_x)))
        # axis(1, at=min(cpu), lab=c(0)) # plot the zero in the beginning

        # draw logarithmic y axis
        axis(2, at=markers_y, lab=markers_y, tck=-.025, cex.axis=0.95, mgp = c(3, 0.5, 0))
    }


}

mtext("Mean percent of booked CPUs used", 1, -.25, outer=TRUE, cex=cex)
mtext("Mean percent of booked RAM used", 2, outer=TRUE, las=0, cex=cex)

if(plot_opt_save_to_file == 1){
	# write to file
	dev.off()
}

