######
# sun angle correlated with tidal height change by site

library(oce)
library(foreach)

cl <- makeCluster(8)
registerDoParallel(cl)

# objects for selecting window width combinations for case studies

cases <- c('ELKVM', 'PDBBY', 'RKBMB', 'SAPDC')
days <- seq(1, 30, by = 1)
its <- expand.grid(days, cases)
its$nms <- paste(its$Var1, its$Var2)
out_ls <- vector('list', nrow(its))
names(out_ls) <- its$nms
for(i in 1:nrow(its)){
  
  # log
  sink('C:/Users/mbeck/Desktop/log.txt')
  cat(i, 'of', nrow(its))
  sink()
  
  # values to iterate
  case <- its$Var2[i]
  daywin <- its$Var1[i]
  
  # file to upload, doesn't matter which
  sel_vec <- paste0(case, '_wtreg_1.RData')
  load(paste0('wtreg/', sel_vec))
  dat <- get(gsub('\\.RData$', '', sel_vec))
  
  # sun angle
  locs <- as.numeric(get_map_meta(case)[, c('Longitude', 'Latitude')])
  utc_time <- as.POSIXlt(dat$DateTimeStamp, tz = 'UTC')
  sun_angle <- sunAngle(utc_time, locs[1], locs[2])$altitude
  
  # tidal change
  tide_change <- dat$Tide
  
  # weights
  dec_time <- dec_fun(dat)
  dec_time <- dec_time[, ncol(dec_time)]
  dat_in <- data.frame(dec_time, hour = dat$hour, Tide = dat$Tide)
  
  cor_out <- foreach(row = 1:nrow(dat)) %dopar% {
    
    ref_in <- dat_in[row, ]
    
    wts <- wt_fun(ref_in, dat_in, wins = list(daywin, 1, 0.6))
    gr_zero <- which(wts > 0)
    
    sun_in <- (sun_angle)[gr_zero]
    tide_in <- (tide_change)[gr_zero]
    
    cor(sun_in, tide_in)
    
  }
 
  cor_out <- unlist(cor_out)
  
  out_ls[[its$nms[i]]] <- cor_out
  
}

angle_tide_its <- out_ls
save(angle_tide_its, file = 'C:/Users/mbeck/Desktop/angle_tide_its1.RData')

tmp <- do.call('cbind', angle_tide_its)

tmp <- data.frame(DateTimeStamp = dat$DateTimeStamp, tmp)
tmp <- melt(tmp, id.var = 'DateTimeStamp')
tmp$daywin <- gsub('^X|\\..*$', '' ,tmp$variable)
tmp$case <- gsub('^X[0-9]*\\.', '', tmp$variable)

pdf('C:/Users/mbeck/Desktop/angle_tide_its1.pdf', height = 9, width = 5)
for(daywin in unique(tmp$daywin)){
  
  p1 <- ggplot(tmp[tmp$daywin %in% daywin, ], 
      aes(x = DateTimeStamp, y = value)) + 
    geom_line() + 
    facet_wrap(~ variable, ncol = 1) + 
    theme_bw() + 
    ylab('Correlation') +
    scale_y_continuous(limits= c(-1, 1)) +
    geom_hline(yintercept = 0, linetype = 'dashed', colour = 'red') +
    ggtitle(paste('Day window', daywin))

  print(p1)
  
}
dev.off()


library(oce)
library(foreach)

cl <- makeCluster(8)
registerDoParallel(cl)

# objects for selecting window width combinations for case studies

cases <- c('ELKVM', 'PDBBY', 'RKBMB', 'SAPDC')
days <- seq(1, 30, by = 1)
its <- expand.grid(days, cases)
its$nms <- paste(its$Var1, its$Var2)
out_ls <- vector('list', nrow(its))
names(out_ls) <- its$nms
for(i in 1:nrow(its)){
  
  # log
  sink('C:/Users/mbeck/Desktop/log.txt')
  cat(i, 'of', nrow(its))
  sink()
  
  # values to iterate
  case <- its$Var2[i]
  daywin <- its$Var1[i]
  
  # file to upload, doesn't matter which
  sel_vec <- paste0(case, '_wtreg_1.RData')
  load(paste0('wtreg/', sel_vec))
  dat <- get(gsub('\\.RData$', '', sel_vec))
  
  # sun angle
  locs <- as.numeric(get_map_meta(case)[, c('Longitude', 'Latitude')])
  utc_time <- as.POSIXlt(dat$DateTimeStamp, tz = 'UTC')
  sun_angle <- sunAngle(utc_time, locs[1], locs[2])$altitude
  
  # tidal change
  tide_change <- dat$Tide
  
  # weights
  dec_time <- dec_fun(dat)
  dec_time <- dec_time[, ncol(dec_time)]
  dat_in <- data.frame(dec_time, hour = dat$hour, Tide = dat$Tide)
  
  cor_out <- foreach(row = 1:nrow(dat)) %dopar% {
    
    ref_in <- dat_in[row, ]
    
    wts <- wt_fun(ref_in, dat_in, wins = list(daywin, 1e6, 1e6))
    gr_zero <- which(wts > 0)
    
    sun_in <- (sun_angle)[gr_zero]
    tide_in <- (tide_change)[gr_zero]
    
    cor(sun_in, tide_in)
    
  }
 
  cor_out <- unlist(cor_out)
  
  out_ls[[its$nms[i]]] <- cor_out
  
}

angle_tide_its <- out_ls
save(angle_tide_its, file = 'C:/Users/mbeck/Desktop/angle_tide_its2.RData')

tmp <- do.call('cbind', angle_tide_its)

tmp <- data.frame(DateTimeStamp = dat$DateTimeStamp, tmp)
tmp <- melt(tmp, id.var = 'DateTimeStamp')
tmp$daywin <- gsub('^X|\\..*$', '' ,tmp$variable)
tmp$case <- gsub('^X[0-9]*\\.', '', tmp$variable)

pdf('C:/Users/mbeck/Desktop/angle_tide_its2.pdf', height = 9, width = 5)
for(daywin in unique(tmp$daywin)){
  
  p1 <- ggplot(tmp[tmp$daywin %in% daywin, ], 
      aes(x = DateTimeStamp, y = value)) + 
    geom_line() + 
    facet_wrap(~ variable, ncol = 1) + 
    theme_bw() + 
    ylab('Correlation') +
    scale_y_continuous(limits= c(-1, 1)) +
    geom_hline(yintercept = 0, linetype = 'dashed', colour = 'red') +
    ggtitle(paste('Day window', daywin))

  print(p1)
  
}
dev.off()


######
# repeat for simulated ts

# grid of simulation conditions to evaluate
tide_cat <- c('Diurnal')
tide_cat <- factor(tide_cat, levels = tide_cat)
bio_rng <- 2 # round(seq(0, 2, length = 3),2)
tide_assoc <- 2 # round(seq(0, 2, length = 3), 2)
err_rng_pro <- 0 # round(seq(0, 2, length = 3), 2)
err_rng_obs <- 0 # round(seq(0, 2, length = 3), 2)

eval_grd <- expand.grid(tide_cat, bio_rng, tide_assoc, err_rng_pro, 
  err_rng_obs)
names(eval_grd) <- c('tide_cat', 'bio_rng', 'tide_assoc', 'err_rng_pro', 
  'err_rng_obs')

dy_wins <- seq(1, 12, by = 1)

# grid of both eval and window combinations
comb_grd <- expand.grid(tide_cat, bio_rng, tide_assoc, err_rng_pro, err_rng_obs,
  dy_wins)
names(comb_grd) <- c(names(eval_grd), 'dy_wins')

# time vector
vec <- c('2014-05-01 00:00:00', '2014-05-31 00:00:00')
vec <- as.POSIXct(vec, format = '%Y-%m-%d %H:%M:%S')
vec <- seq(vec[1], vec[2], by = 60*30)

# setup parallel 
cl <- makeCluster(8)
registerDoParallel(cl)

# iterate through evaluation grid to create sim series
# saved as separate files in /prdnrm/
strt <- Sys.time()
out_ls <- foreach(row = 1:nrow(comb_grd)) %dopar% {
 
  # progress
  sink('C:/Users/mbeck/Desktop/log.txt')
  cat(row, ' of ', nrow(comb_grd), '\n')
  sink()
    
  # eval grid to evaluate
  to_eval <- comb_grd[row, ]
  
  # create simulated time series of DO, tide, etc.
  DO_sim <- with(to_eval, 
    ts_create(
      vec, 
      do.amp = bio_rng, 
      tide_cat = as.character(tide_cat), 
      tide_assoc = tide_assoc,
      err_rng_obs = err_rng_obs,
      err_rng_pro = err_rng_pro
      )  
    )
  
  # corr windows
  win_in <- list(comb_grd[row, 'dy_wins'], 1e6, 1e6)

  cors <- NULL
  for(i in 1:nrow(DO_sim)){
    
    ref_in <- DO_sim[i, ]
    
    wts <- wt_fun(ref_in, DO_sim, wins = win_in)
    gr_zero <- which(wts > 0)
    
    sun_in <- (DO_sim$DO_bio)[gr_zero]
    tide_in <- (DO_sim$Tide)[gr_zero]
    
    cors <- c(cors, cor(sun_in, tide_in))
    
  }
  
  cors
  
}
names(out_ls) <- paste('dy_wins', comb_grd$dy_wins) 

tmp <- do.call('cbind', out_ls)

tmp <- data.frame(DateTimeStamp = vec, tmp)
tmp <- melt(tmp, id.var = 'DateTimeStamp')
tmp$daywin <- gsub('^dy_wins\\.', '' ,tmp$variable)

pdf('C:/Users/mbeck/Desktop/sim_corrs.pdf', height = 4, width = 8)
for(daywin in unique(tmp$daywin)){
  
  p1 <- ggplot(tmp[tmp$daywin %in% daywin, ], 
      aes(x = DateTimeStamp, y = value)) + 
    geom_line() + 
    theme_bw() + 
    ylab('Correlation') +
    scale_y_continuous(limits= c(-1, 1)) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    ggtitle(paste('Day window', daywin))

  print(p1)
  
}
dev.off()

