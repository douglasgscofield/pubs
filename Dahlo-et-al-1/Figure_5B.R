# $ module load R/3.2.3

# example run:
# Rscript Figure_5B.R '2016-01-01' '2017-01-01'


# rm(list=ls())
# library('RMySQL')
library("RSQLite")
library('stringr')
library('RColorBrewer')
library('MASS')

# for devel use
dateStart = "2016-01-01"
dateEnd = "2017-01-01"

# get arguments
args = args<-commandArgs(TRUE)
dateStart = args[1]
dateEnd = args[2]


# connect to the db
jobs_db = dbConnect(dbDriver("SQLite"),"../data/efficiency.sqlite")

	
# print(dateStart)
# print(dateEnd)


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
mem = with(jobs, tapply(mem_peak *counts, proj_id, FUN = sum))
data = data.frame(pid = names(mem), mem = mem, cpu = with(jobs, tapply(cpu_mean*cores *counts, proj_id, FUN = sum)), corecounts = with(jobs, tapply(corecounts, proj_id, FUN = sum)), counts = with(jobs, tapply(counts, proj_id, FUN = sum)), stringsAsFactors = F, row.names=NULL)


# calculate the average again
data$mem = data$mem/data$counts
data$cpu = data$cpu/data$counts




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


# find all the a projects
lv = is.na(str_match(jobs$proj_id, "^plat\\d+$"))
jobs$type[!lv] = groups[3]

# find all the b projects
lv = is.na(str_match(jobs$proj_id, "^bio\\d+$"))
jobs$type[!lv] = groups[2]

# find all the uppmax projects, ie not a or b projs
lv = is.na(str_match(jobs$proj_id, "^gen\\d+$"))
jobs$type[!lv] = groups[1]



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




    ### numbers used in the paper under Results and Discussion - Resource usage
		
		
    ### numbers used in the paper under Results and Discussion - Efficiency

    # prepare data
    jobs$core_usage = jobs$cpu_mean * jobs$cores

    # percentage using 1 core or less
    dim(jobs[jobs$type=='UPPNEX' & jobs$core_usage <= 100,])[1] / dim(jobs[jobs$type=='UPPNEX',])[1]
	dim(jobs[jobs$type=='UPPMAX' & jobs$core_usage <= 100,])[1] / dim(jobs[jobs$type=='UPPMAX',])[1]



    # cpu quantiles
    # 1+ core
    quantile(jobs$core_usage[jobs$type=='UPPNEX' & jobs$core_usage > 100])
    quantile(jobs$core_usage[jobs$type=='UPPMAX' & jobs$core_usage > 100])
    quantile(jobs$core_usage[jobs$type=='PLATFORM' & jobs$core_usage > 100])

    # <1 core
    quantile(jobs$core_usage[jobs$type=='UPPNEX' & jobs$core_usage <= 100])
    quantile(jobs$core_usage[jobs$type=='UPPMAX' & jobs$core_usage <= 100])
    quantile(jobs$core_usage[jobs$type=='PLATFORM' & jobs$core_usage <= 100])


    # percentage using more than 16 cores
    dim(jobs[jobs$type=='UPPNEX' & jobs$core_usage > 1600,])[1] / dim(jobs[jobs$type=='UPPNEX',])[1]
	dim(jobs[jobs$type=='UPPMAX' & jobs$core_usage > 1600,])[1] / dim(jobs[jobs$type=='UPPMAX',])[1]




    # memory quantiles
    # 1+ core
    quantile(jobs$mem_peak[jobs$type=='UPPNEX' & jobs$core_usage > 100])
    quantile(jobs$mem_peak[jobs$type=='UPPMAX' & jobs$core_usage > 100])
    quantile(jobs$mem_peak[jobs$type=='PLATFORM' & jobs$core_usage > 100])

    # <1 core
    quantile(jobs$mem_peak[jobs$type=='UPPNEX' & jobs$core_usage <= 100])
    quantile(jobs$mem_peak[jobs$type=='UPPMAX' & jobs$core_usage <= 100])
    quantile(jobs$mem_peak[jobs$type=='PLATFORM' & jobs$core_usage <= 100])



    # number of jobs with an average usage above 16 cores
    dim(jobs[jobs$core_usage>=1600 & jobs$type=='UPPNEX',])[1]
    dim(jobs[jobs$core_usage>=1600 & jobs$type=='UPPNEX',])[1]/dim(jobs[jobs$type=='UPPNEX',])[1]
    dim(jobs[jobs$core_usage>=1600 & jobs$type=='UPPMAX',])[1]
    dim(jobs[jobs$core_usage>=1600 & jobs$type=='UPPMAX',])[1]/dim(jobs[jobs$type=='UPPMAX',])[1]

    # quote of the percentages above
    (dim(jobs[jobs$core_usage>=1600 & jobs$type=='UPPMAX',])[1]/dim(jobs[jobs$type=='UPPMAX',])[1])/(dim(jobs[jobs$core_usage>=1600 & jobs$type=='UPPNEX',])[1]/dim(jobs[jobs$type=='UPPNEX',])[1])



    # same as above, but with core hours instead
    sum(jobs$corecounts[jobs$core_usage>=1600 & jobs$type=='UPPNEX'])
    sum(jobs$corecounts[jobs$core_usage>=1600 & jobs$type=='UPPNEX'])/sum(jobs$corecounts[jobs$type=='UPPNEX'])
    sum(jobs$corecounts[jobs$core_usage>=1600 & jobs$type=='UPPMAX'])
    sum(jobs$corecounts[jobs$core_usage>=1600 & jobs$type=='UPPMAX'])/sum(jobs$corecounts[jobs$type=='UPPMAX'])

    # quote of the percentages above
    (sum(jobs$corecounts[jobs$core_usage>=1600 & jobs$type=='UPPMAX'])/sum(jobs$corecounts[jobs$type=='UPPMAX']))/(sum(jobs$corecounts[jobs$core_usage>=1600 & jobs$type=='UPPNEX'])/sum(jobs$corecounts[jobs$type=='UPPNEX']))


    # not correct to take mean/median directly, should factor in counts
    # mean(jobs$mem_peak[jobs$type=='UPPNEX' & jobs$counts>1])
    # median(jobs$mem_peak[jobs$type=='UPPNEX' & jobs$counts>1])
    # mean(jobs$mem_peak[jobs$type=='UPPMAX' & jobs$counts>1])
    # median(jobs$mem_peak[jobs$type=='UPPMAX' & jobs$counts>1])

    # counts factored in, but only for calculating mean. Median would involve doing a indexes like 20 lines down
    # sum(jobs$cpu_mean[jobs$type=='UPPMAX' & jobs$counts>1] * jobs$cores[jobs$type=='UPPMAX' & jobs$counts>1] * jobs$counts[jobs$type=='UPPMAX' & jobs$counts>1]) / sum(jobs$counts[jobs$type=='UPPMAX' & jobs$counts>1])
    # sum(jobs$cpu_mean[jobs$type=='UPPNEX' & jobs$counts>1] * jobs$cores[jobs$type=='UPPNEX' & jobs$counts>1] * jobs$counts[jobs$type=='UPPNEX' & jobs$counts>1]) / sum(jobs$counts[jobs$type=='UPPNEX' & jobs$counts>1])

    # like above, but not updated to include counts as above
    # median(jobs$cpu_mean[jobs$type=='UPPNEX' & jobs$counts>1])
    # mean(jobs$cpu_mean[jobs$type=='UPPMAX' & jobs$counts>1])
    # median(jobs$cpu_mean[jobs$type=='UPPMAX' & jobs$counts>1])

}


weighted_data = jobs
weighted_data$cpu = weighted_data$cpu_mean * weighted_data$cores
weighted_data$corecounts = weighted_data$cores * weighted_data$counts

cpu = weighted_data$cpu
mem = weighted_data$mem_peak


# convert memory usage to megabytes
mem = mem*1024


### adjust data for logarithmic axis
# set ~zero mem usage to 100 megabyte to compress the plot
mem[mem<128] = 128

# set the ~zero cpu usage to 50% of a core to compress the plot
cpu[cpu<50] = 50

# convert cpu usage to core usage
cpu = cpu/100

# logaritmize the values
log_base = 2
cpu_log = log(cpu, base=log_base)
mem_log = log(mem, base=log_base)


### create the surface plot

# bins
nbins <- 1 # production
# nbins <- 300 # devel

# KernSmooth kde with abkde estimated bandwidth
surf <- lapply(split(data.frame(cpu_log[1:1000], mem_log[1:1000]), weighted_data$type[1:1000]),
    function(x){ 
   
            bandwidth=c(1.5,1.5)
            kde2d(x[,1], x[,2], n=nbins, h=bandwidth)
        })

# rescale the heatmap intensities
orgsurf = surf
# surf = orgsurf
for(i in 1:length(surf)){
	temp = data.frame(surf[[i]])
	surf[[i]]$z = data.matrix(temp[,3:dim(temp)[2]]/100)



}

# make the axis go from 0 - 101 to make it even
ax = seq(min(cpu_log) , max(cpu_log), length.out=nbins+1)
ay = seq(min(mem_log) , max(mem_log), length.out=nbins+1)



# for both uppmax and uppnex
names = c('Non-NGS projects', 'NGS projects', 'Platform')
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
    types = c(4:8) # De novo, RNA, Methylation, WGS, Metagenomics
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
plot_opt_ncol = 100 # the number of colors the heat map should have
plot_opt_image = 0 # plot the heat map?
plot_opt_scale_image = 0 # scale the heat man intensity after group size?
plot_opt_points = 1 # plot the data points?
plot_opt_subplot_titles = 1 # plot title for each subplot?
plot_opt_contour = 0 # plot the contours in the heat map?
plot_opt_grid = 1 # plot a grid in the plot?
plot_opt_axis = 1 # plot axes in the plot?
plot_opt_save_to_file = 1 # save to a file?# plot uppmax and uppnex
plot_opt_file_name = 'Figure_5b_resource_usage.png'


if(plot_opt_save_to_file == 1){
	# save to file

	xy = c(4.33, 2.55)
	r = 1200
    png(plot_opt_file_name, xy[1]*r, xy[2]*r, res=r)
    cex = 0.7
    par(mfrow=c(1,length(plot_order)), las=1, mar=c(2, 2, 2, .5), oma=c(1, 1.25, 0, .25), cex=cex)
}else{
    cex=2.2
    par(mfrow=c(1,length(plot_order)), las=1, mar=c(2, 3, 3, .5), oma=c(3, 2, 1, 1), cex=cex)
}


# set the plot properties, 2 sub-plots, margins around in outer and inner area




for(i in 1:length(plot_order)){

    # pick out the index of surf that contains the project type we are interested in
    surf_idx = which(names(surf)==plot_order[i])
    g_idx = which(groups==plot_order[i])

    
    	
	# generate colors
	col='#00000002'

	# 'increase the exposure'
	# col = c(col, rep(col[length(col)], 5000))
	
	if(plot_opt_scale_image == 1){
    	# scale the colorRampPalette after the size of category. 
    	col = col[1: ceiling(plot_opt_ncol * count_sum[surf_idx] / max(count_sum)) ]
    }
    	
    	

    if(plot_opt_image == 1){
        # plot the surface
        image(surf[[surf_idx]], col=col , useRaster=TRUE, xlab="", ylab="", axes=FALSE, ylim=c(floor(log2(min(mem))), ceiling(log2(max(mem))) ), xlim=c(floor(log2(min(cpu))), ceiling(log2(max(cpu))) ) )
        
    }else{
    	plot(c(0,0), col='#FFFFFF00', xlab="", ylab="", axes=FALSE, ylim=c(floor(log2(min(mem))), ceiling(log2(max(mem))) ), xlim=c(floor(log2(min(cpu))), ceiling(log2(max(cpu))) ) )
    }



    if(plot_opt_points == 1){
        # plot data points
        subsample = as.integer(runif(20000,1,dim(weighted_data)[1]))
        subsample = 1:dim(weighted_data)[1]
        points(cpu_log[weighted_data$type == groups[g_idx]][subsample], mem_log[weighted_data$type == groups[g_idx]][subsample], cex=.32, pch=16, col=col)
        # text(cpu_log[weighted_data$type == groups[i]], mem_log[weighted_data$type == groups[i]], weighted_data$pid[weighted_data$type == groups[i]], cex=0.5)
    }

    if(plot_opt_subplot_titles == 1){
        # add title to subplot
        title(names[g_idx])
    }


    if(plot_opt_contour == 1){
        # draw a contour
        contour(surf[[surf_idx]], add=T, lwd=.3, drawlabels=FALSE, nlevels=10)
    }


    # create axis marks
    markers_x = (1:ceiling(max(cpu)))
    markers_x = (1:2^ceiling(log(max(cpu))/log(2)))
    markers_x = markers_x[log(markers_x, base=log_base) %% 1 == 0]
    markers_y = c(1:ceiling(log2(max(mem))))

    if(plot_opt_grid == 1){
        # print grid lines
        abline( v = c(0, log(markers_x, base=log_base)), lty = 3, col = brewer.pal(9, "Greys")[4], lwd=.3 )
        abline( h = markers_y, lty = 3, col = brewer.pal(9, "Greys")[4], lwd=.3)
    }


    if(plot_opt_axis == 1){
        # draw logarithmic x axis (subtract 0.005 to get the last mark in the plot)
        # axis(1, lab=FALSE, tck=-.025)
        axis(1, at=c(min(cpu_log),log(markers_x, base=log_base)), lab=FALSE, tck=-.025) # plot the zero in the beginning

        text(c(min(cpu_log),log(markers_x, base=log_base)), par("usr")[3] - 0.65, srt = 45, adj = 1, labels = c(0,markers_x), xpd = TRUE, cex=.85)


        # draw logarithmic y axis
        # axis(2, at=markers_y, lab=(log_base^markers_y)/1024)
        axis(2, at=markers_y[7:length(markers_y)], lab=c(0,((log_base^markers_y)/1024)[8:length(markers_y)]), tck=-.025, cex.axis=0.85, mgp = c(3, 0.5, 0) ) # hmm, 2048 gb not visible.. 
        # mgp=c(1, 3, 1)
    }


}

mtext("Mean CPU usage", 1, -.25, outer=TRUE, cex=cex)
mtext("Max GiB of RAM used", 2, 0.25, outer=TRUE, las=0, cex=cex)

if(plot_opt_save_to_file == 1){
	# write to file
	dev.off()
}


