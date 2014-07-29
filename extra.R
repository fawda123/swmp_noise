# ######
# # figure representing data in tables 1/2
# 
# load('mod_perf.RData')
# 
# # calculate mean, sd for prd_cor and prd_err
# prd_perf <- llply(names(mod_perf)[1:8],
#   .fun = function(x) {
#     ddply(
#       mod_perf,
#       .variables = x,
#       .fun = function(y) {
#         c(
#           mean(y$prd_cor),
#           sd(y$prd_cor),
#           mean(y$prd_err),
#           sd(y$prd_err)
#         )
#     })}
#   )
# # rename first column for combining, then combine
# prd_perf <- lapply(prd_perf, 
#   function(x) { 
#     Parameter <- names(x)[1]
#     names(x)[1] <- 'Level'
#     x[,1] <- factor(x[,1])
#     x[,1] <- c('Lo', 'Med', 'High')
#     Parameter <- rep(Parameter, 3)
#     x <- data.frame(Parameter, x)
#     x
#     }
#   )
# prd_perf <- do.call('rbind', prd_perf)
# names(prd_perf) <- c('Parameter', 'Level', 'Mean cor', 'Sd cor', 'Mean err', 'Sd err')
# 
# mean_prd_perf <- melt(prd_perf, id.var = c('Parameter', 'Level'), 
#   measure.var = c('Mean cor', 'Mean err'))
# mean_prd_perf$variable <- gsub('Mean ', '', mean_prd_perf$variable)
# sd_prd_perf <- melt(prd_perf, id.var = c('Parameter', 'Level'), 
#   measure.var = c('Sd cor', 'Sd err'))
# sd_prd_perf$variable <- gsub('Sd ', '', sd_prd_perf$variable)
# 
# to_plo <- merge(mean_prd_perf, sd_prd_perf, 
#   by = c('variable', 'Parameter', 'Level'))
# names(to_plo)[grep('value', names(to_plo))] <- c('Mean', 'Sd')
# labs <- c('Lo', 'Med', 'High')
# to_plo$Level <- factor(to_plo$Level, labels = labs, levels = labs)
# to_plo$Eval <- 'Predicted'
# 
# ####
# 
# # calculate mean, sd for dtd_cor and dtd_err
# dtd_perf <- llply(names(mod_perf)[1:8],
#   .fun = function(x) {
#     ddply(
#       mod_perf,
#       .variables = x,
#       .fun = function(y) {
#         c(
#           mean(y$dtd_cor),
#           sd(y$dtd_cor),
#           mean(y$dtd_err),
#           sd(y$dtd_err)
#         )
#     })}
#   )
# # rename first column for combining, then combine
# dtd_perf <- lapply(dtd_perf, 
#   function(x) { 
#     Parameter <- names(x)[1]
#     names(x)[1] <- 'Level'
#     x[,1] <- c('Lo', 'Med', 'High')
#     x$Parameter <- rep(Parameter, 3)
#     x
#     }
#   )
# dtd_perf <- do.call('rbind', dtd_perf)
# names(dtd_perf) <- c('Level', 'Mean cor', 'Sd cor', 'Mean err', 'Sd err', 'Parameter')
# 
# mean_dtd_perf <- melt(dtd_perf, id.var = c('Parameter', 'Level'), 
#   measure.var = c('Mean cor', 'Mean err'))
# mean_dtd_perf$variable <- gsub('Mean ', '', mean_dtd_perf$variable)
# sd_dtd_perf <- melt(dtd_perf, id.var = c('Parameter', 'Level'), 
#   measure.var = c('Sd cor', 'Sd err'))
# sd_dtd_perf$variable <- gsub('Sd ', '', sd_dtd_perf$variable)
# 
# to_plo2 <- merge(mean_dtd_perf, sd_dtd_perf, 
#   by = c('variable', 'Parameter', 'Level'))
# names(to_plo2)[grep('value', names(to_plo2))] <- c('Mean', 'Sd')
# labs <- c('Lo', 'Med', 'High')
# to_plo2$Level <- factor(to_plo2$Level, labels = labs, levels = labs)
# to_plo2$Eval <- 'Detided'
# 
# to_plo_all <- rbind(to_plo, to_plo2)
# 
# p1 <- ggplot(to_plo_all, aes(x = Parameter, y = Mean, fill = Level)) +  
#   geom_bar(stat = 'identity', position = position_dodge(.9)) + 
#   geom_errorbar(position = position_dodge(.9), width=.25, 
#     aes(ymin = Mean - Sd, ymax = Mean + Sd)) +
#   facet_wrap(Eval ~ variable, scales = 'free_y') +
#   theme_bw() +
#   theme(legend.title = element_blank()) +
#   scale_fill_grey() 
# 
# ######
# # examples of how detiding works, use for shiny
# 
# # combination grid, used for systematic sims
# load('comb_grd.RData')
# 
# # load predicted/normalized data
# load('prdnrm.RData')
# 
# # get results from prdnrm corresponding to a given window value in which_wins
# # rownames in which_wins are used to iterate through prdnrm
# names(comb_grd)[names(comb_grd) %in% 'dec_time'] <- 'dec_win'
# 
# sel_vec <- with(comb_grd,
#   bio_rng %in% 2 &
#   tide_assoc %in% 0 &
#   err_rng_obs %in% 0 &
#   err_rng_pro %in% 0
#   )
# 
# # sub of comb_grd 
# comb_tmp <- comb_grd[sel_vec, ]
# comb_tmp$L1 <- rownames(comb_tmp)
# # sub of prdnrm list, name list with rownames of comb_grd for merge
# prd_tmp <- prdnrm[sel_vec]
# names(prd_tmp) <- rownames(comb_tmp)
# 
# # melt predictions and merge with comb_tmp
# prd_tmp <- melt(prd_tmp, id.var = names(prd_tmp[[1]]))
# prd_tmp <- merge(prd_tmp, comb_tmp, by = 'L1', all.x = T)
# 
# # create residual vector
# prd_tmp$DO_res <- with(prd_tmp, DO_obs - DO_prd)
# # dobule windows since theyre halves
# prd_tmp$dec_win <- as.character(2 * prd_tmp$dec_win)
# 
# levs <- c('Tide', 'DO_bio', 'DO_adv', 'DO_obs', 'DO_prd',
#   'DO_nrm', 'DO_dtd')
# to.plo <- melt(prd_tmp, id.var = c('Day', 'dec_win', 'tide_cat'),
#   measure.var = levs
#   )
# to.plo$variable <- factor(to.plo$variable, levels = levs)
# levs <- as.character(c(2, 10, 20, 30, 40))
# to.plo$dec_win <- factor(to.plo$dec_win, labels = levs, levels = levs)
# 
# to.plos <- split(to.plo, to.plo$tide_cat)
# 
# ylab<-expression(paste(O[2], ' (mg ',L^-1,')'))
# for(p in 1:length(to.plos)){
# pls <- ggplot(to.plos[[p]], aes(x = Day, y = value, 
#     group = dec_win, colour = dec_win)) +
#   geom_line(size = 1) +
#   facet_wrap(~ variable, scales = 'free_y', ncol = 1) + 
#   theme_bw() +
#   ylab(ylab) + 
#   scale_colour_manual('Window', values = brewer.pal(5, 'Spectral')) +
#   theme(legend.position = 'top') +
#   ggtitle(unique(to.plos[[p]]$tide_cat))
# assign(paste0('p', p), pls)
# }
# grid.arrange(p1, p2, p3, ncol = 3)
# 
# # facet_wrap_labeller(p, labels = c(
# #   expression(italic(DO [obs])),
# #   expression(italic(DO [tid])),
# #   expression(italic(DO [mtd])),
# #   expression(italic(
# #     paste(DO [res], '=', DO [obs] - DO [tid]))),
# #   expression(italic(
# #     paste(DO [dtd], '=', DO [mtd] + DO [res])))
# #   ))
# 
# ######
# # moving window correlations??

load('met_ls_inst.RData')

# subs <- format(met_ls_inst[[1]]$DateTimeStamp, '%m') %in% '07'

cl <- makeCluster(4)
registerDoParallel(cl)

win_cors <- ldply(met_ls_inst, 
  .parallel = T, 
  .fun = function(x) {
    
    corrs <- matrix(NA_real_, nrow = nrow(x), ncol = 4)
    for(i in 1:nrow(x)) {
      
      cat(i, '\t')
  
      ref_in <- x[i, ]
  
      wts <- wt_fun(ref_in, x, win = 2)
      sel_vec <- wts > 0
  
      DO_obs <- try(with(x[sel_vec, ], 
          cor(DO_obs, Tide, use = 'complete.obs')      
          ))
      if('try-error' %in% class(DO_obs)) DO_obs <- NA
        
      DO_dtd <- try(with(x[sel_vec, ], 
          cor(DO_dtd, Tide, use = 'complete.obs')      
          ))
      if('try-error' %in% class(DO_dtd)) DO_dtd <- NA
      
      DOF_obs <- try(with(x[sel_vec, ], 
          cor(DOF_obs, dTide, use = 'complete.obs')      
          ))
      if('try-error' %in% class(DOF_obs)) DOF_obs <- NA
      
      DOF_dtd <- try(with(x[sel_vec, ], 
          cor(DOF_dtd, dTide, use = 'complete.obs')      
          ))
      if('try-error' %in% class(DOF_dtd)) DOF_dtd <- NA
      
      
      corrs[i, ] <- c(DO_obs, DO_dtd, DOF_obs, DOF_dtd)  
      
      }
    
    out <- data.frame(x[, c('DateTimeStamp', 'Tide', 'dTide')], corrs)
    names(out)[grep('X', names(out))] <- c('DO_obs', 'DO_dtd', 'DOF_obs', 'DOF_dtd')
    
    return(out)

    }
  
  )

to_plo <- melt(win_cors, 
  measure.var = c('DO_obs', 'DO_dtd', 'DOF_obs', 'DOF_dtd')
  )
to_plo$wtreg <- gsub('DO_|DOF_', '', to_plo$variable )
to_plo$variable <- gsub('_obs|_dtd', '', to_plo$variable)

ggplot(to_plo, aes(x = DateTimeStamp, y = value, colour = wtreg)) + 
  geom_line() + 
  theme_bw() +
  facet_grid(.id ~ variable)

######
# table of DO/metab correlations before after, detiding
# note that tide in met_ls is daily average of hourly tidal change

# metab and inst flux data
load('met_ls.RData')
load('met_ls_inst.RData')
load('case_grds.RData')

# go through each site for DO cors, use metab list for metab cors
case_regs <- list.files(getwd(), '_wtreg_[0-9]*.RData')
cor_res <- alply(matrix(case_regs),
  1, 
  .fun = function(x){
  
    # load wtreg data
    load(x)
    nm <- gsub('.RData', '', x)
    dat_in <- get(nm)
      
    # DO obs v tide
    do_obs <- with(dat_in, 
      cor.test(DO_obs, Tide)
      )
    
    # DO dtd v tide
    do_dtd <- with(dat_in, 
      cor.test(DO_dtd, Tide)
      )
    
    # inst flux data
    flux_in <- met_ls_inst[[x]]
    
    # DOF with dtide
    flux_obs <- with(flux_in, 
      cor.test(DOF_obs, dTide)
      )
    
    # DOF_dtd with dtide
    flux_dtd <- with(flux_in, 
      cor.test(DOF_dtd, dTide)
      ) 
    
    # get tidal range for metabolic day/night periods from flux_in
    # for correlation with daily integrated metab
    tide_rngs <- ddply(flux_in, 
      .variables = c('met.date'),
      .fun = function(x){
#         sunrise <- suppressWarnings(diff(range(x[x$variable %in% 'sunrise', 'Tide'])))
#         sunset <- suppressWarnings(diff(range(x[x$variable %in% 'sunset', 'Tide'])))
#         if(sunrise == 'Inf') sunrise <- NA
#         if(sunset == 'Inf') sunset <- NA
#         daytot <- diff(range(x$Tide))
        sunrise <- mean(diff(x[x$variable %in% 'sunrise', 'Tide'], na.rm = T))
        sunset <- mean(diff(x[x$variable %in% 'sunset', 'Tide'], na.rm = T))
        if(sunrise == 'Inf') sunrise <- NA
        if(sunset == 'Inf') sunset <- NA
        daytot <- mean(diff(x$Tide, na.rm = T))
        
        c(daytot, sunrise, sunset)
        }
      )
    names(tide_rngs) <- c('Date','daytot', 'sunrise', 'sunset')
    
    # get metab data from list
    dat_in <- met_ls[[x]]
    dat_in <- merge(dat_in, tide_rngs, by = 'Date', all.x = T)
    
    # as list for all metab correlations
    # Pg values correlated with tidal range during sunlight hours
    # Rt values correlated with tidal range during night hours
    # NEM values correlated with metabolic daily tidal range
    met_cor <- list(
      
      Pg_obs = with(dat_in, 
        cor.test(Pg, sunrise)
        ),
    
      Rt_obs = with(dat_in, 
        cor.test(Rt, sunset)
        ),
    
      NEM_obs = with(dat_in, 
        cor.test(NEM, daytot)
        ),
    
      Pg_dtd = with(dat_in, 
        cor.test(Pg_dtd, sunrise)
        ),
    
      Rt_dtd = with(dat_in, 
        cor.test(Rt_dtd, sunset)
        ),
    
      NEM_dtd = with(dat_in, 
        cor.test(NEM_dtd, daytot)
        )
      
      )
    
    # DO and metab corrs combined
    all_ls <- c(do_obs = list(do_obs), do_dtd = list(do_dtd),
      flux_obs = list(flux_obs), flux_dtd = list(flux_dtd), met_cor)
    
    # convert the stats for each wtreg to data frame
    res_sum <- ldply(all_ls, 
      function(x) with(x, c(estimate, p.value))
      )
    names(res_sum) <- c('var', 'cor', 'pval')

    res_sum
    
  })
names(cor_res) <- case_regs

# melt and make separate columns for site and window comb value
cor_res <- melt(cor_res, id.var = names(cor_res[[1]]))  
cor_res$site <- gsub('_wtreg_[0-9]*.RData', '', cor_res$L1)
cor_res$wins <- as.numeric(gsub('^.*_wtreg_|.RData', '', cor_res$L1))

# merge with case_grds
case_grds$wins <- as.numeric(row.names(case_grds))
cor_res <- merge(cor_res, case_grds, by = 'wins', all.x = T)

# create columns for variable (DO, flux, etc.) and sub variable (obs, dtd)
cor_res$sub_var <- gsub('^.*_', '', cor_res$var)
cor_res$var <- gsub('_.*$', '', cor_res$var)

save(cor_res, file = 'cor_res.RData')

to_plo <- cor_res
to_plo$group_var <- paste(to_plo$Tide, to_plo$sub_var)
to_plo_obs <- to_plo[to_plo$sub_var %in% 'obs', ]
p1 <- ggplot(to_plo[to_plo$sub_var %in% 'dtd',], 
    aes(x = factor(Day), y = cor, colour = Tide, group = group_var)) +
  geom_line() + 
  geom_line(data = to_plo_obs, 
    aes(x = factor(Day), y = cor, group = group_var), 
    colour = 'black', size = 1) +
  geom_point(aes(pch = sub_var)) +
  facet_grid(var ~ site) +
  ylim(c(-1, 1)) +
  theme_bw() 

  