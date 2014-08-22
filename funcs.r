#functions for NEM processing of NERRS data
#created Dec. 2013 by M. Beck, adapted from 'spam_NEM_fun.r' and M. Murrell

#funcion that splits dataset into 24hr days based on sunrise
#merge with original data
met.day.fun<-function(dat.in, stat.in, 
  meta.path = 'M:/wq_models/SWMP/sampling_stations.csv'
  ){

  require(StreamMetabolism)  #for sunrise.set function
  
  if(!exists('dat.meta')) 
    dat.meta<-read.csv(meta.path,header=T)
  
  stat.meta<-toupper(paste0(stat.in,'WQ'))
  stat.meta<-dat.meta[grep(stat.meta,toupper(dat.meta$Station.Code)),]
  
  # all times are standard - no DST!
  gmt.tab<-data.frame(
    gmt.off=c(-4,-5,-6,-8,-9),
    tz=c('America/Virgin', 'America/Jamaica', 'America/Regina',
      'Pacific/Pitcairn', 'Pacific/Gambier'),
    stringsAsFactors=F
    )
  
  #get sunrise/sunset times using sunrise.set function from StreamMetabolism
  lat<-stat.meta$Latitude
  long<--1*stat.meta$Longitude
  GMT.Offset<-stat.meta$GMT.Offset
  tz<-gmt.tab[gmt.tab$gmt.off==GMT.Offset,'tz']
  start.day<-format(dat.in$DateTimeStamp[which.min(dat.in$DateTimeStamp)]-(60*60*24),format='%Y/%m/%d')
  tot.days<-1+length(unique(as.Date(dat.in$DateTimeStamp)))
  
  #ss.dat is matrix of sunrise/set times for each days  within period of obs
  ss.dat<-suppressWarnings(sunrise.set(lat,long,start.day,tz,tot.days))
  
  #remove duplicates, sometimes sunrise.set screws up
  ss.dat<-ss.dat[!duplicated(strftime(ss.dat[,1],format='%Y-%m_%d')),]
  ss.dat<-data.frame(
    ss.dat,
    met.date=as.Date(ss.dat$sunrise,tz=tz)
    )
  ss.dat<-melt(ss.dat,id.vars='met.date')
  if(!"POSIXct" %in% class(ss.dat$value))
    ss.dat$value<-as.POSIXct(ss.dat$value, origin='1970-01-01',tz=tz)
  ss.dat<-ss.dat[order(ss.dat$value),]
  ss.dat$day.hrs<-unlist(lapply(
    split(ss.dat,ss.dat$met.date),
    function(x) rep(as.numeric(x[2,'value']-x[1,'value']),2) 
    ))
  
  #matches is vector of row numbers indicating starting value that each
  #unique DateTimeStamp is within in ss.dat
  #output is meteorological day matches appended to dat.in
  matches<-findInterval(dat.in$DateTimeStamp,ss.dat$value)
  data.frame(dat.in,ss.dat[matches,])
      
  }

######
#gets station name, first site then station, separated by comma
#input is five character string for site uppe case, e.g., 'ACEBB'
stat.name.fun<-function(stat){
  
  if(!exists('dat.meta')) 
    dat.meta<-read.csv('M:/wq_models/SWMP/sampling_stations.csv',header=T,
      strip.white=T)
  
  out<-grep(tolower(stat),dat.meta$Station.Code)
  out<-unique(dat.meta[out,names(dat.meta) %in% c('Station.Name','Reserve.Name')])
  out<-paste(out[,'Reserve.Name'],out[,'Station.Name'],sep=', ')
  
  return(out)
  
  }

######
#creates simulated DO curve from sine wave
#input is vector of POSIX date objects at 30 min interals
#'do.mean' is used to change intercept
#'do.amp' changes amplitude
#'freq' changes period, this value is set to the number of observations within a period
DO.bio.fun<-function(time.in,sunrise='06:30:00',do.mean=8,do.amp=1,freq=48){
 
  # time.step per hour
  time.step <- 60/unique(diff(time.in))
  
  #time diff between strt of vector and sunrise
  sunrise <- paste(as.Date(min(time.in)), sunrise)
  sunrise <- as.POSIXct(sunrise, tz = attr(time.in, 'tzone'))
  cor.phs <- time.step * as.numeric(difftime(sunrise, min(time.in)))
  cor.phs <- 2 * pi * (1/freq) * cor.phs
  cor.phs <- cor.phs - pi
	
  # get input values to cos func
	in.vals <- seq(0,length(time.in), length = length(time.in))
  in.vals <- in.vals * 2 * pi * (1/freq)

	#amplitude change
	amp<-do.amp
	
	#intercept
	int<-do.mean
  
  # DO signal
	DO<-int+amp*cos(in.vals - cor.phs)
	
	return(DO)
	
	}

######
# creates simulated tidal signal as additive combination of sine waves
# resulting tide is scaled from 0 to 1 and added to 'scl.up'
# uses DO.bio.fun above
# 'time.in' is posix time vector for simulating tides
# 'waves.in' is list of each tidal component with two elements, period and amp
# 'scl.up' value to add to additive tidal vector
# 'all' is logical if each component should be returned
tide.fun <- function(time.in, waves.in, scl.up = 4, all = F){
  
  require(scales)
  
  # simulate each component from DO.bio.fun
  tide <- sapply(
    waves.in, 
    function(x){
      x[1] <- 2*x[1] # for 30 min obs
      DO.bio.fun(time.in, do.mean = 0, 
        do.amp = x[[2]], freq = x[[1]])
      }
    ) 
  
  # return all components if T
  if(all) return(tide)
  
  # add each component, center/scale at zero, 
  # change range to 1m, shift up to some arbitrary mean height
  tide <- scl.up + scales::rescale(scale(rowSums(tide)))
  
  # combine with data frame and add dtide/dt
  dTide <- c(diff(tide)[1], diff(tide))        

  # output
  out <- data.frame(DateTimeStamp = time.in, Tide = tide, dTide = dTide)

  return(out)
  
  }

######
# calculate number of anomolous estimates from nem output
# 'nem.in' is output from 'nem.fun'
# 'pgvar' and 'rtvar' are character strings of metab values to use
anoms.fun <- function(nem.in, pgvar = 'Pg', rtvar = 'Rt') {
  
  Pg <- nem.in[, pgvar]
  Pg <- sum(Pg <= 0,na.rm = T)/length(na.omit(Pg))
  Rt <- nem.in[, rtvar]
  Rt <- sum(Rt >= 0,na.rm = T)/length(na.omit(Rt))
  
  out <- data.frame(Pg,Rt)
  names(out) <- c(pgvar, rtvar)
  
  return(out)
  
  }

######
# calculate number of anomolous estimates from flux input
# 'inst.in' is output from 'inst.flux.fun'
# 'stat.in' is text string of station to eval, depends on 'met.day.fun'
# 'flxvar' is character of flux variable to estime
anoms.inst.fun <- function(inst.in, stat.in, flxvar = 'DOF') {
  
  # get metabolic day info
  to_eval <- met.day.fun(inst.in, stat.in)
  
  # separate records into day/night
  flx_day <- to_eval[to_eval$variable %in% 'sunrise', flxvar]
  flx_nght <- to_eval[to_eval$variable %in% 'sunset', flxvar]
  
  # get anoms by day (negative flux), night (positive flux)
  anom_day <- sum(flx_day <= 0, na.rm = T)/length(flx_day)
  anom_nght <- sum(flx_nght >= 0, na.rm = T)/length(flx_nght)
  
  # output
  out <- data.frame(anom_day, anom_nght)
  
  return(out)
  
  }

######
#r.squared function
#created for data from weighted regression, but works for all models
#residuals are observed - predicted, but taken from model objects (see 'epc_mods.R')
rsq.fun<-function(resid,obs){
  
  require(Metrics)
  
  ssr<-sum(resid^2,na.rm=T)
  sst<-sum(se(obs,mean(obs,na.rm=T)),na.rm=T)
  
  return(1 - (ssr/sst))
  
}

######
#variant of rmse fun in Metrics package but handles na values
#resid is obs - predicted
rmse.fun<-function(resid){
  
  out<-sqrt(mean(resid^2,na.rm=T))
    
  return(out) 
  
  }

######
# found on SO
# creates summary data frame from a larger data frame
# 'data' is input df
# 'measurevar' is variable to summarize
# 'groupvars' is chr vector of variables to condition summary
# similar to aggregate function
summarySE <- function(data=NULL, measurevar, groupvars=NULL, narm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    require(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=narm) {
        if (narm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=narm),
          mean = mean   (xx[[col]], na.rm=narm),
          sd   = sd     (xx[[col]], na.rm=narm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

######
# function for plotting repeating polygons in ggplot
# 'flag.in' is vector indicating binomial variable for polys
# 'flag.in' can be factor or numeric (w/ two values)
# 'dat' is data.frame with 'flag.in'
# 'for_leg' is used to manually change color for polygons, will add legend
# output is geom object
poly.fun<-function(flag.in,dat, for_leg = F){

  require(reshape2)
  require(ggplot2)
  
  #for flag bias
  if(class(flag.in) == 'numeric'){ 
    neg.dates<-with(
      dat,
      DateTimeStamp[which(flag.in < 0)]
      )
    tz<-attr(neg.dates,'tzone')
    diffs<-c(0,diff(as.numeric(neg.dates)))
    strt.ind<-c(1,which(diffs>1800))
    end.ind<-c(strt.ind-1,length(flag.in))
    comb<-paste(neg.dates[strt.ind],neg.dates[end.ind],sep="\t")
    
    if(grepl('NA',comb[length(comb)])) 
      comb[length(comb)]<-gsub('NA',neg.dates[length(neg.dates)],comb[length(comb)])
  
    comb<-do.call('rbind',strsplit(comb,'\t'))
    comb<-cbind(comb,comb[,2],comb[,1])
  
    x.vals<-suppressMessages(melt(sapply(1:nrow(comb), 
      function(x) comb[x,],simplify=F))$value)
    x.vals<-as.POSIXct(as.character(x.vals),tz,
      format='%Y-%m-%d %H:%M:%S')
    y.vals<-rep(c(-1000,-1000,1000,1000),length=length(x.vals)) 
    Antagonistic<-rep(1:(length(x.vals)/4),each=4)
    polys<-data.frame(x.vals,y.vals,grp=Antagonistic)
    
    }

  #for sunset/rise
  if(class(flag.in) == 'factor'){
    
    plo.dates<-unique(dat[,c('solar','value')])
    
    if(plo.dates$solar[1] == 'sunset') 
      plo.dates<-plo.dates[-1,]
    if(plo.dates$solar[nrow(plo.dates)] == 'sunrise')
      plo.dates<-rbind(
        plo.dates,
        data.frame(solar='sunset',value=max(dat$DateTimeStamp))
        )
    
    plo.dates$inds<-rep(1:(nrow(plo.dates)/2),each=2)
    tz<-attr(plo.dates$value,'tzone')
    plo.dates$value <- as.character(plo.dates$value)
    plo.dates<-dcast(plo.dates,inds~solar,value.var='value')
    plo.dates<-with(plo.dates,
      data.frame(sunrise,sunset,sunset,sunrise)
      )
   
    x.vals <- sapply(1:nrow(plo.dates), 
           function(x) plo.dates[x,],simplify=F)
    x.vals<-suppressMessages(
      melt(x.vals, measure.vars = names(x.vals[[1]]))$value
      )
    x.vals<-as.POSIXct(x.vals, tz, origin = '1970-01-01')
    y.vals<-rep(c(-1000,-1000,1000,1000),nrow(plo.dates))
    Day<-as.character(trunc(x.vals,'days'))
    polys<-data.frame(x.vals,y.vals,grp=Day)

    }
  
  if(for_leg){
    out<-geom_polygon(data=polys,aes(x.vals,y.vals,group=grp, fill = 'grp'), 
      alpha=0.6)
  } else {
   out<-geom_polygon(data=polys,aes(x.vals,y.vals,group=grp), fill = 'orange',
      alpha=0.6)
  }
    
  
  return(out)
  
  }

######
# create dec time using jday on 24 hour scale
# 'dat_in' is data frame input with time vecot as 'var_nm'
# output is same data frame including 'dec_time' column
# note that this is different from the one used for cases
dec_fun <- function(dat_in, var_nm = 'DateTimeStamp'){

  # get time vec from dat_in
  posix_in  <- dat_in[, var_nm]
  
  #dec time on 24 hr scale
  jday <- as.numeric(format(posix_in, '%j')) 
  hour <- as.numeric(format(posix_in, '%H'))
  minu<- as.numeric(format(posix_in, '%M'))
  
  # add separate vecs on same scale
  dec_time <- jday + (hour/24) + (minu/60/24)
  
  # convert hour back to 24 hour scale, add minutes
  hour <- hour + minu/60
  
  out <- data.frame(dat_in, jday, hour, dec_time)
  return(out)
  
  }
######
# this is a repliate of filled.contour that removes category borders in the legend
# http://stackoverflow.com/questions/8068366/removing-lines-within-filled-contour-legend
filled.contour.hack <- function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, 
    length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
    ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
    levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
    col = color.palette(length(levels) - 1), plot.title, plot.axes, 
    key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
    axes = TRUE, frame.plot = axes, ...) 
{
    if (missing(z)) {
        if (!missing(x)) {
            if (is.list(x)) {
                z <- x$z
                y <- x$y
                x <- x$x
            }
            else {
                z <- x
                x <- seq.int(0, 1, length.out = nrow(z))
            }
        }
        else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
        y <- x$y
        x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
        stop("increasing 'x' and 'y' values expected")
    mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    on.exit(par(par.orig))
    w <- (3 + mar.orig[2L]) * par("csi") * 2.54
    layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
    par(las = las)
    mar <- mar.orig
    mar[4L] <- mar[2L]
    mar[2L] <- 1
    par(mar = mar)
    plot.new()
    plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
        yaxs = "i")
    rect(0, levels[-length(levels)], 1, levels[-1L], col = col, border = NA)
    if (missing(key.axes)) {
        if (axes) 
            axis(4)
    }
    else key.axes
    box()
    if (!missing(key.title)) 
        key.title
    mar <- mar.orig
    mar[4L] <- 1
    par(mar = mar)
    plot.new()
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    .filled.contour(x, y, z, levels, col)
    if (missing(plot.axes)) {
        if (axes) {
            title(main = "", xlab = "", ylab = "")
            Axis(x, side = 1)
            Axis(y, side = 2)
        }
    }
    else plot.axes
    if (frame.plot) 
        box()
    if (missing(plot.title)) 
        title(...)
    else plot.title
    invisible()
  }

######
# creates data for DO simulation, uses functions above
# 'time_in' is POSIX vector of continuous time series to estimate, currently as 30 minute obs
# 'do.amp' is amplitude of biological DO signal
# 'tide_cat' is either a list of tidal components, 
# a character string as 'Diurnal', 'Semidiurnal', or 'Mixed Semidiurnal', 
# or a data.frame of tidal values
# 'tide_assoc' is value +- range of DO_adv when converting Tide (i.e., amp)
# 'err_rng_obs' is sd of obs uncertainty for DO signal
# 'err_rng_pro' is bounded range +/- or proc unc
# 'seeded' is logical indicating if a random seed is used for obs/pro error
# output is data drame with DateTimeStamp, jday, hour, dec_time,
# sunrise, Tide, DO_bio, DO_adv, DO_tid, DO_noi, DO_obs
ts_create <- function(time_in, do.amp, tide_cat, tide_assoc, err_rng_obs, 
  err_rng_pro, seeded = T){

  require(scales)
  
  # start data.frame for plotting
  DO_sim <- data.frame(
    DateTimeStamp = time_in
    )
  
  # add decimal time on julian day
  DO_sim <- dec_fun(DO_sim)

  # add column for dec_time starting from one (for plotting)
  DO_sim$Day <- with(DO_sim, (dec_time + 1) - jday[1])
  
  # add vector for sunrise, -1 is daylight, 1 is nighttime
  DO_sim$sunrise <- 1
  day_vec <- with(DO_sim, 24*(dec_time - jday))
  DO_sim$sunrise[day_vec > 6 & day_vec <= 18] <- -1
  
  # move dec_time to correspond to meteorological day
  # note this is lazy, only works if sunrise is 630
  DO_sim$dec_time  <- DO_sim$dec_time - 6.5/24
  
  # get tide if character
  if(grepl('character|factor', class(tide_cat))){
    if(tide_cat == 'Semidiurnal')
      waves.in <- list(c(12.42, 1)) # M2, principal lunar semidiurnal
    if(tide_cat == 'Diurnal') 
      waves.in <- list(c(25.82, 1)) #O1, principal lunar diurnal
    if(tide_cat == 'Mixed Semidiurnal')
      waves.in <- list(c(12.42, 1), c(25.82, 1)) 
    
    tide <- tide.fun(vec, waves.in)
    }

  # get tide if list
  if(class(tide_cat) == 'list')
    tide <- tide.fun(vec, waves.in)
  
  # get tide if supplied as data.frame, first col is posix, second is tide
  if(class(tide_cat) == 'data.frame'){
    tide <- tide_cat
    tide$Tide <- scales::rescale(tide$Tide, to = c(4, 5))
    }
  
  # combine with data frame and add dtide/dt
  DO_sim$Tide <- tide$Tide
  
  ##
  # create biological DO
  
  # obs and proc uncertainty
  if(seeded) set.seed(1234)
  DO_sim$e_obs <- rnorm(nrow(DO_sim), 0, err_rng_obs)
  if(seeded) set.seed(4321)
  DO_sim$e_pro <- scales::rescale(cumsum(rnorm(nrow(DO_sim), 0, err_rng_pro)), 
    to = c(-1 * err_rng_pro, err_rng_pro))
  DO_sim$DO_unc <- with(DO_sim, e_obs + e_pro)

  # add uncertainty to clean DO signal
  DO_sim$DO_die <- DO.bio.fun(time_in, do.amp = do.amp)
  DO_sim$DO_bio <- pmax(0, with(DO_sim, DO_die + DO_unc))

  # add tidal effect
  DO_sim$DO_adv <- scales::rescale(DO_sim$Tide, to = c(-1*tide_assoc, tide_assoc))
  DO_sim$DO_obs <- with(DO_sim, DO_bio + DO_adv)

  # floor at zero
  DO_sim$DO_obs <- with(DO_sim, pmax(0, DO_obs))
  
  # 
  
  return(DO_sim)
  
  }

######
# using expressions in facet_wrap labels
# http://stackoverflow.com/questions/19282897/how-to-add-expressions-to-labels-in-facet-wrap
facet_wrap_labeller <- function(gg.plot,labels=NULL) {
  #works with R 3.0.1 and ggplot2 0.9.3.1
  require(gridExtra)

  g <- ggplotGrob(gg.plot)
  gg <- g$grobs      
  strips <- grep("strip_t", names(gg))

  for(ii in seq_along(labels))  {
    modgrob <- getGrob(gg[[strips[ii]]], "strip.text", 
                       grep=TRUE, global=TRUE)
    gg[[strips[ii]]]$children[[modgrob$name]] <- editGrob(modgrob,label=labels[ii])
  }

  g$grobs <- gg
  class(g) = c("arrange", "ggplot",class(g)) 
  g
}

######
# gets station location data for ggmap plotting
# 'val_in' is character string of first chrs of reserve, or five chrs for a site
# if only reserve is given, all sites are returned
# input is not case-sensitive
# returns appropriate row(s) from sampling_stations.csv
get_map_meta <- function(val_in){
  
  #removing trailing whit space in chr strings
  trim.trailing <- function(x) sub('^\\s+|\\s+$', '', x)

  # station metadata
  stats<-read.csv('M:/wq_models/SWMP/sampling_stations.csv',header=T)
  stats<-stats[grep('Active*',stats$Status),]
  
  # subset by reserve
  stats <- stats[grep(paste0(tolower(val_in),'.*wq'), stats[,'Station.Code']),]
  stats$Longitude <- -1*stats$Longitude
  stats$Station <- trim.trailing(as.character(stats$Station.Name))
  stats$Station.Code <- toupper(substr(trim.trailing(as.character(stats$Station.Code)),1,5))
  
  return(stats)
  
  }

######
# some stupid function for formatting p-values for tables
p_form<-function(x, tex_out=F){
  if(x<0.005) out<-'< 0.005'
  else{
    if(x<0.05&x>=0.005) out<-as.character(format(round(x,3),nsmall=3))
    else out<-'ns'
    }
  if(tex_out == F) out
  else{
    if(out=='< 0.005') '$p<$ 0.005'
    else{
      if(out=='ns') '$p=$ ns'
      else paste('$p=$',out,sep=' ')
      }
    }
  }

######
# formatting of values in S expressions
form_fun <- function(x, rnd_val = 2, dig_val = 2, nsm_val = 2) {
  format(round(x, rnd_val), digits = dig_val, nsmall = nsm_val)
  }

######
# removing trailing whit space in chr strings
# 'x' is character string
trim.trailing <- function(x) sub('^\\s+|\\s+$', '', x)

######
# gets station location data for ggmap plotting
# 'val_in' is character string of first chrs of reserve, or five chrs for a site
# if only reserve is given, all sites are returned
# input is not case-sensitive
# returns appropriate row(s) from sampling_stations.csv
get_map_meta <- function(val_in){
  
  # station metadata
  stats<-read.csv('M:/wq_models/SWMP/sampling_stations.csv',header=T)
  stats<-stats[grep('Active*',stats$Status),]
  
  # subset by reserve
  stats <- stats[grep(paste0(tolower(val_in),'.*wq', collapse = '|'), stats[,'Station.Code']),]
  stats$Longitude <- -1*stats$Longitude
  stats$Station <- trim.trailing(as.character(stats$Station.Name))
  stats$Station.Code <- toupper(substr(trim.trailing(as.character(stats$Station.Code)),1,5))
  
  return(stats)
  
  }

######
# creates significance stars from vector of p values
# 'dat_in' is vector of p-values
star_fun <- function(dat_in){
  
  if(!'numeric' %in% class(dat_in)) 
    stop('Input data must be numeric vector')
  
  symnum(dat_in,
    cutpoints=c(.001,0.01,.05),
    corr=T,
    symbols=c("***","**","*"," "),
    abbr=T,diag=F) 
  
  }

######
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

######
# get predicted, normalized values not using interp grid, tide as predictor
# 'dat_in' is raw data used to create 'grd_in' and used to get predictions
# 'DO_obs' is string indicating name of col for observed DO values from 'dat_in'
# output is data frame same as 'dat_in' but includes predicted and norm columns
wtreg_fun <- function(dat_in, DO_obs = 'DO_obs', wins = list(4, 12, NULL),
  parallel = F, progress = F){

  # get mean tidal height from empirical data
  mean_tide <- mean(dat_in$Tide)

  #for counter
  strt <- Sys.time()
  
  out <- ddply(dat_in, 
    .variable = 'DateTimeStamp',
    .parallel = parallel, 
    .fun = function(row){
      
      # row for prediction
      ref_in <- row
      ref_in <- ref_in[rep(1, 2),]
      ref_in$Tide <- c(unique(ref_in$Tide), mean_tide)
      
      # progress
      if(progress){
        prog <- which(row$DateTimeStamp == dat_in$DateTimeStamp)
        sink('log.txt')
        cat('Log entry time', as.character(Sys.time()), '\n')
        cat(prog, ' of ', nrow(dat_in), '\n')
        print(Sys.time() - strt)
        sink()
        }
      
      # get wts
      ref_wts <- wt_fun(ref_in, dat_in, wins = wins, slice = T, 
        subs_only = T, wt_vars = c('dec_time', 'hour', 'Tide'))
  
      #OLS wtd model
      out <- lapply(1:length(ref_wts),
        function(x){
          
          # subset data for weights > 0
          dat_proc <- dat_in[as.numeric(names(ref_wts[[x]])),]
          
          # if no DO values after subset, return NA
          # or if observed DO for the row is NA, return NA
          if(sum(is.na(dat_proc$DO_obs)) == nrow(dat_proc)|
              any(is.na((ref_in$DO_obs)))){
            
            DO_pred <- NA
            beta <- NA
            Tide <- ref_in$Tide[x]
            
            } else {
            
              # subset weigths > 0, rescale weights average
              ref_wts <- ref_wts[[x]]/mean(ref_wts[[x]])
            
              # get model
              mod_md <- lm(
                DO_obs ~ dec_time + Tide, # + sin(2*pi*dec_time) + cos(2*pi*dec_time),
                weights = ref_wts,
                data = dat_proc
                )
            
              # get prediction from model
              Tide <- ref_in$Tide[x]
              DO_pred <- predict(
                mod_md, 
                newdata = data.frame(dec_time = ref_in$dec_time[x], Tide = Tide)
                )
            
              # get beta from model
              beta <- mod_md$coefficients['Tide']
            
            }
          
          # output
          DO_pred
          
          }
        
        )

      out <- unlist(out)
      names(out) <- c('DO_prd', 'DO_nrm')
      out
      
      })
  
  out$DateTimeStamp <- NULL
  out <- cbind(dat_in, out)

  return(out)
  
  }

######
#function for getting regression weights
# note that this subsets the input data frame for faster wt selection
# subset is by limiting window for product of weights (dec_time)
# subsetted weights are recombined to equal vector of length = nrow(dat_in)
#'wt_vars' is name of three variables to weight
#'ref_in' is row of dat.in that is used as reference
#'dat_in' is data to get weights from
#'wins' are the windows for the three wt.vars, values represent halves
#'all' will return all weights, rather than the product of all three
#'slice' is logical for subsetting 'dat_in' for faster wt selection
#'subs_only' is logical for returning only wt vectors that are non-zero
wt_fun <- function(ref_in, dat_in,
  wt_vars = c('dec_time', 'hour', 'Tide'),
  wins = list(4, 12, NULL),
  all = F, 
  slice = T, 
  subs_only = F){
  
  # sanity check
  if(sum(wt_vars %in% names(dat_in)) != length(wt_vars))
    stop('Weighting variables must be named in "dat_in"')
  
  # windows for each of three variables
  wins_1<-wins[[1]]
  wins_2<-wins[[2]]
  wins_3<-wins[[3]]
  
  # default window width for third variable is half its range
  if(is.null(wins[[3]])) wins_3 <- diff(range(dat_in[, wt_vars[3]]))/2
  
  # weighting tri-cube function
  # mirror extends weighting function if vector repeats, e.g. monthly
  # 'dat_cal' is observation for weight assignment
  # 'ref' is reference observation for fitting the model
  # 'win' is window width from above (divided by two)
  # 'mirr' is logical indicating if distance accounts for repeating variables (e.g., month)
  # 'scl_val' is range for the ref vector of obs, used to get correct distance for mirrored obs
  wt_fun_sub <- function(dat_cal, ref, win, mirr = F, scl_val = 1){
    
    # dist_val is distance of value from the ref
    dist_val <- sapply(ref, function(x) abs(dat_cal - x))
    
    # repeat if distance is checked on non-continuous number line
    if(mirr){
      
        dist_val <- pmin(
          sapply(ref, function(x)
            abs(x + scl_val - dat_cal)),
          sapply(ref, function(x) abs(dat_cal + scl_val - x)),
          dist_val
          )
      
      }
    
    # get wts within window, otherwise zero
    win_out <- dist_val > win
    dist_val <- (1 - (dist_val/win)^3)^3
    dist_val[win_out] <- 0
      
    return(dist_val)
      
    }

  #reference (starting) data
  ref_1 <- as.numeric(ref_in[, wt_vars[1]])
  ref_2 <- as.numeric(ref_in[, wt_vars[2]])
  ref_3 <- as.numeric(ref_in[, wt_vars[3]])

  ##
  # subset 'dat_in' by max window size for faster calc
  # this is repeated if min number of wts > 0 is not met
  # subset vector is all T if not using subset
  dec_rng <- range(dat_in$dec_time)
  ref_time <- unique(ref_in$dec_time)
  dec_sub <- with(dat_in, 
    dec_time > 
      ref_time - wins_1 * 5 & dec_time < ref_time + wins_1 * 5
    )
  if(!slice) dec_sub <- rep(T, length = nrow(dat_in))
  dat_sub <- dat_in[dec_sub, ]

  ##
  # weights for each observation in relation to reference
  # see comments for 'wt_fun_sub' for 'scl_val' argument
  
  # jday
  wts_1 <- wt_fun_sub(as.numeric(dat_sub[, wt_vars[1]]), 
    ref = ref_1, win = wins_1, mirr = F) 
  # hour
  wts_2 <- wt_fun_sub(as.numeric(dat_sub[, wt_vars[2]]), 
    ref = ref_2, win = wins_2, mirr = T, scl_val = 24)
  # tide
  wts_3 <- wt_fun_sub(as.numeric(dat_sub[, wt_vars[3]]), 
    ref = ref_3, win = wins_3, mirr = F)
  # all as product 
  out <- sapply(1:nrow(ref_in), function(x) wts_1[, x] * wts_2[, x] * wts_3[, x])
  
  gr_zero <- colSums(out > 0)
  #cat('   Number of weights greater than zero =',gr.zero,'\n')
  
  # extend window widths of weight vector is less than 100
  while(any(gr_zero < 100)){
    
    # increase window size by 10%
    wins_1 <- 1.1 * wins_1
    wins_2 <- 1.1 * wins_2
    wins_3 <- 1.1 * wins_3 
    
    # subset again
    dec_sub <- with(dat_in, 
      dec_time > ref_time - wins_1 * 5 & dec_time < ref_time + wins_1 * 5
      )
    if(!slice) dec_sub <- rep(T, length = nrow(dat_in))
    dat_sub <- dat_in[dec_sub, ]
    
    #weights for each observation in relation to reference
    wts_1 <- wt_fun_sub(as.numeric(dat_sub[, wt_vars[1]]), 
      ref = ref_1, win = wins_1, mirr = F)
    wts_2 <- wt_fun_sub(as.numeric(dat_sub[, wt_vars[2]]), 
      ref = ref_2, win = wins_2, mirr = T, scl_val = 24)
    wts_3 <- wt_fun_sub(as.numeric(dat_sub[, wt_vars[3]]), 
      ref = ref_3, win = wins_3, mirr = F)
    
    out <- sapply(1:nrow(ref_in), 
      function(x) wts_1[, x] * wts_2[, x] * wts_3[, x])
    
    gr_zero <- colSums(out > 0)
    
    }
  
  if(subs_only){
    
    nms <- which(dec_sub)
    out <- alply(out, 2, function(x) {
    
      to_sel <- x > 0
      tmp <- x[to_sel]
      names(tmp) <- which(dec_sub)[to_sel]
      tmp
    
      })
    
    return(out)
    
    }
  
  # extend weight vectors to length of dat_in
  empty_mat <- matrix(0, ncol = nrow(ref_in), nrow = nrow(dat_in))
  empty_fill <- function(wts_in) {
    out <- empty_mat
    out[dec_sub,] <- wts_in
    out
    }
  wts_1 <- empty_fill(wts_1)
  wts_2 <- empty_fill(wts_2)
  wts_3 <- empty_fill(wts_3)  
  out <- empty_fill(out)

  #return all weights if T
  if(all){
    out <- data.frame(dat_in$DateTimeStamp, 
      wts_1, wts_2, wts_3, out)
    names(out) <- c('DateTimeStamp', wt_vars, 'final')
    return(out)
    }
  
  #final weights are product of all three
  out
  
  }

######
# get legend from an existing ggplot object
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

######
# returns default ggplot colors
ggplotColours <- function(n=6, h=c(0, 360) +15){
  if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}