
# example run to crunch the numbers, with recalculate = 1 below:
# Rscript Figure_7.R 2 7 '2010-07-01' '2017-01-01'
# Rscript Figure_7.R 3 7 '2010-07-01' '2017-01-01'

# example run to make the plots, with recalculate = 0 below:
# Rscript Figure_7.R 2 7 '2010-07-01' '2017-01-01'


# rm(list=ls())
library("RSQLite")
library("RColorBrewer")

# connect to the db
jobs_db = dbConnect(dbDriver("SQLite"),"../data/efficiency.sqlite")

# should the number be recalculated again? Only if the algorithm has changed or the db has been updated.
recalculate = 0

# get arguments
args = args<-commandArgs(TRUE)
type = as.integer(args[1])
timeperiod = as.integer(args[2])
dateStart = args[3]
dateEnd = args[4]





# the part below collects the data and calculates everything. It saves the workspace at the end of the if statement, so unless you change anything in the number crunching you should just skip the whole if statement and load the workspace instead, as it does on lint ~250.
if(recalculate){

	library("Rcpp")



	# for(type in 1:4){

	if(type == 1){
		# get platform corehour usage per week
		results <- dbSendQuery(jobs_db, 
		paste("select * from jobs where counts > 1 and proj_id like 'plat%' and date_finished>='",dateStart,"' and date_finished<'",dateEnd,"';", sep='')) # avoid jobs that are less than 10 minutes long
		jobs = fetch(results,n=-1) # get all results
		dbClearResult(results) # reset db
		plot_filename = "eot_plat"
	}



	if(type == 2){
		# get uppnex corehour usage per week
		results <- dbSendQuery(jobs_db, 
		paste("select * from jobs where counts > 1 and (proj_id like 'bio%') and date_finished>='",dateStart,"' and date_finished<'",dateEnd,"';", sep='')) # avoid jobs that are less than 10 minutes long
		jobs = fetch(results,n=-1) # get all results
		dbClearResult(results) # reset db
		plot_filename = "eot_ngs"
	}


	if(type == 3){
		# get uppmax corehour usage per week
		results <- dbSendQuery(jobs_db, 
		paste("select * from jobs where counts > 1 and (proj_id like 'gen%') and date_finished>='",dateStart,"' and date_finished<'",dateEnd,"';", sep='')) # avoid jobs that are less than 10 minutes long
		jobs = fetch(results,n=-1) # get all results
		dbClearResult(results) # reset db
		plot_filename = "eot_gen"
	}


	### for readability

	# convert dates to Date
	finished_date = as.Date(jobs$date_finished)

	# calculate how many days the job could strech
	runtime_in_days = as.integer(jobs$counts/288) # efficiency is sampled once every 5 min -> 288 times per day

	# Set the start date of the job
	start_date = finished_date - runtime_in_days

	# divide the counts per day the job ran
	# counts_per_day = jobs$counts / (runtime_in_days+1)

	# save to the data.frame
	jobs$date_finished = finished_date
	jobs$date_started = start_date

	jobs$week_finished = strftime(as.POSIXlt(jobs$date_finished),format="%W")

	jobs$runtime_in_days = runtime_in_days+1
	#jobs$counts_per_day = counts_per_day


	# adjust memory usage above the limit to avoid >100% 'efficiency'
	jobs$mem_peak = pmin(jobs$mem_peak, jobs$mem_limit)


	# first, calculate the actual efficiency volume (cores x counts x efficiency)
	jobs$mem_efficiency = 100*jobs$mem_peak/jobs$mem_limit *jobs$counts*jobs$cores 
	jobs$cpu_efficiency = jobs$cpu_mean *jobs$counts*jobs$cores

	# the jobs efficiency is the max of either cpu or mem
	jobs$job_efficiency = pmax(jobs$mem_efficiency, jobs$cpu_efficiency)

	# then calculate the the max volume of the job (the volume it would have had if it was 100% efficient all the time)
	jobs$max_efficiency = 100.0 *jobs$counts*jobs$cores 


	# create a sequence of days between the days (skipping the last day since it will have no jobs running on it, just ending)
	days = as.data.frame(seq(min(jobs$date_started), max(jobs$date_finished)-1, "days"))
	colnames(days) = c('date')


	Sys.time()
	# init
	#percentiles = list()
	means = c()
	medians = c()
	p25 = c()
	p75 = c()

	means_efficiency = c()
	medians_efficiency = c()
	p25_efficiency = c()
	p75_efficiency = c()

	# for each week in the time line
	for(day in days$date[seq(1,length(days$date), timeperiod)]){

		day = as.Date(day, origin="1970-01-01")
		print(day)


		# get all jobs covering this day
		period_jobs = jobs[jobs$date_finished>=day & jobs$date_started<(day+timeperiod),]

		# skip to the next week if there are no jobs in this one
		if(dim(period_jobs)[1] == 0){
			#percentiles = c(percentiles, NA)
			means = c(means, NA)
			medians = c(medians, NA)
			p25 = c(p25, NA)
			p75 = c(p75, NA)

			means_efficiency = c(means_efficiency, NA)
			medians_efficiency = c(medians_efficiency, NA)
			p25_efficiency = c(p25_efficiency, NA)
			p75_efficiency = c(p75_efficiency, NA)
			next
		}


		# get how many days of the jobs overlap this week
		pre_overlap = as.numeric(day - period_jobs$date_started)
		pre_overlap[pre_overlap<0] = 0
		post_overlap = period_jobs$date_finished - (day+timeperiod-1)
		post_overlap[post_overlap<0] = 0
		period_jobs$overlap = as.numeric((period_jobs$runtime_in_days - post_overlap - pre_overlap)/period_jobs$runtime_in_days)



		# calculate how many percent each job makes out of the total this week
		period_jobs$percentage = period_jobs$counts*period_jobs$cores*period_jobs$overlap / sum(period_jobs$counts*period_jobs$cores*period_jobs$overlap)

		# preallocate size
		percentile = rep(NA, sum(ceiling(period_jobs$counts*period_jobs$cores*period_jobs$overlap)))
		# generate percentile entries
		counter = 1
		for(i in 1:dim(period_jobs)[1]){
			current_length = ceiling(period_jobs$counts[i]*period_jobs$cores[i]*period_jobs$overlap[i])
			percentile[counter:(counter + current_length)] = 100*period_jobs$job_efficiency[i]/period_jobs$max_efficiency[i]

			counter = counter + current_length

		}

		# save the values
		means = c(means, mean(percentile))
		medians = c(medians, median(percentile))
		p25 = c(p25, quantile(percentile)[2])
		p75 = c(p75, quantile(percentile)[4])

		means_efficiency = c(means_efficiency, mean(period_jobs$job_efficiency/period_jobs$max_efficiency))
		medians_efficiency = c(medians_efficiency, median(period_jobs$job_efficiency/period_jobs$max_efficiency))
		p25_efficiency = c(p25_efficiency, quantile(period_jobs$job_efficiency/period_jobs$max_efficiency)[2])
		p75_efficiency = c(p75_efficiency, quantile(period_jobs$job_efficiency/period_jobs$max_efficiency)[4])




	}
	Sys.time()


	data = data.frame(date=days$date[seq(1,length(days$date), timeperiod)], mean=means, median=medians, p25=p25, p75=p75, mean_efficiency=means_efficiency*100, median_efficiency=medians_efficiency*100, p25_efficiency=p25_efficiency*100, p75_efficiency=p75_efficiency*100, col = "#EF3B2C")
	data = data[!is.na(data$mean),]

	# remove date limits
	rm(dateStart)
	rm(dateEnd)
	recalculate = 0

	save.image(paste("percentiles.", type, ".", timeperiod, ".RData", sep=''))


	stop(paste("Saved workspace: percentiles.", type, ".", timeperiod, ".RData", sep=''))



}




# get a rolling mean of the values
rollmean <- function(x,window_size) {
	window_size = ceiling(window_size/2)
    result = rep(0,length(x))
    for(i in 1:length(result)){
        result[i] = mean(x[max(1, i-window_size):min(length(x), i+window_size)])
    }
    return(result)
}


# get a rolling median of the values
rollmedian <- function(x,window_size) {
	window_size = ceiling(window_size/2)
    result = rep(0,length(x))
    for(i in 1:length(result)){
        result[i] = median(x[max(1, i-window_size):min(length(x), i+window_size)])
    }
    return(result)
}




for(type in c(2,3)){

	load(paste("percentiles.", type, ".", timeperiod, ".RData", sep=''))
	print(paste("NOTE: Possible date span from saved data: ",min(data$date),' - ',max(data$date),sep=''))

	# adjust date range of data
	data = data[data$date<=dateEnd & data$date>=dateStart,]

	### color the dots depending on system

	# set plot options based on which groups it is
	if(plot_filename == 'eot_plat'){
		main="Platform projects"

		# milou
		data$col = brewer.pal(9,name='Blues')[7]

	}else if(plot_filename == 'eot_ngs'){

		# milou
		data$col = brewer.pal(9,name='BuGn')[7]

		main="NGS projects"
		color = 'BuGn'

	}else if(plot_filename == 'eot_gen'){
		main="Non-NGS projects"
		data$col = brewer.pal(9,name='YlOrBr')[6]

	}


	### plot the data
	
	if(type==2){
		pdf(file='Figure_7_efficiency_over_time.pdf', 4.33, 2.35)
	}else if(type==1){
		pdf(file='Figure_7_efficiency_over_time.pdf', 2.165, 2.35)
	}

	cex = 0.5

	if(type==2){
		par(mfrow=c(1,2), las=1, xaxs="i", yaxs="i", cex.lab=cex, mar=c(.75, 1.1, 1.2, .5), oma=c(1.3, 1, .5, .5))

	}else if(type==1){
		par(las=1, xaxs="i", yaxs="i", cex.lab=cex, mar=c(.75, 1.1, 1.2, .7), oma=c(1.3, 1, .5, .5))
	}

	loess_span = 0.65

	
	# plot with median values
	scatter.smooth(x=data$date, y=data$median, col=data$col, ylim=c(0,100), xlim=c(min(data$date), max(max(data$date), as.Date("2017-01-01"))), xaxt = "n", xlab='', ylab="", pch=16, axes=FALSE, span=loess_span, lpars=list(lwd=3), bty="L", cex=0.5)

	grid(nx=NA, ny=NULL, col="#00000050")
	title(main=main, line=.4, cex.main=0.83)

	at = c(seq(as.Date(paste(strftime(min(data$date), format='%Y'), '-01-01', sep='')), as.Date(paste(strftime(max(data$date), format='%Y'), '-01-01', sep='')), by='year'), as.Date(paste(as.integer(strftime(max(data$date), format='%Y'))+1, '-01-01', sep='')))
	labels = c(strftime(seq(as.Date(paste(strftime(min(data$date), format='%Y'), '-01-01', sep='')), as.Date(paste(strftime(max(data$date), format='%Y'), '-01-01', sep='')), by='year'), format='%Y'), as.integer(strftime(max(data$date), format='%Y'))+1)
	axis(1, at=at, labels=labels, tck=-.025, cex.axis=0.63, padj=-2.5)

	# if it's uppnex projects being plotted
	if(type==2){
		# add vertical line to show where efficiency reporting started
		arrows(x0=as.Date("2014-03-31"), y0=15, x1=as.Date("2014-03-31"), y1=0, length=0.06)
		text(as.Date("2014-03-31"), 10, labels=c("Efficiency feedback\n to users began"), adj=-.05, cex=0.6)


	}

	axis(2, cex.axis=0.63, tck=-.025, labels=rep("",6), at=seq(0,100,by=20))
	axis(2, cex.axis=0.63, tick=FALSE, at=seq(0,100,by=20), line=-.65)
	box()

	if(type==2 | type==1){
		mtext("Efficiency (%)", 2, outer=TRUE, las=0, cex=0.7, ) # padj=-6
		mtext("Date", 1, outer=TRUE, las=0, cex=0.7) # , padj=5.3, adj=0
	}

	
	if(type==3){
		dev.off()
	}else if(type==1){
		dev.off()
	}

}

