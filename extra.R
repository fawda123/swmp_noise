
######
# repeat but get monthly correlations?
# metab and inst flux data
load('met_ls.RData')
# load('met_ls_inst.RData')
load('case_grds.RData')

# go through each site for DO cors, use metab list for metab cors
case_regs <- list.files('wtreg', '_wtreg_[0-9]*.RData')
cor_res <- alply(matrix(case_regs),
  1, 
  .progress = 'tk',
  .fun = function(x){
  
    # load wtreg data
    load(paste0('wtreg/', x))
    nm <- gsub('.RData', '', x)
    dat_in <- get(nm)
      
    # add month column
    dat_in$month <- strftime(dat_in$Date, '%m')
    
    # DO obs v tide
    do_obs <- ddply(
      dat_in, 
      .variable = c('month'), 
      .fun = function(x) with(x, cor.test(DO_obs, Tide)$estimate)
    )
    
    # DO dtd v tide
    do_dtd <- ddply(
      dat_in, 
      .variable = c('month'), 
      .fun = function(x) with(x, cor.test(DO_nrm, Tide)$estimate)
    )
    
    # get tidal range for metabolic day/night periods from flux_in
    # for correlation with daily integrated metab
    tide_rngs <- ddply(dat_in, 
      .variables = c('met.date'),
      .fun = function(x){
        
        # mean tidal derivative for day hours
        sunrise <- mean(diff(x[x$variable %in% 'sunrise', 'Tide'], 
          na.rm = T))
        
        # mean tidal derivative for night hours
        sunset <- mean(diff(x[x$variable %in% 'sunset', 'Tide'], 
          na.rm = T))
        if(sunrise == 'Inf') sunrise <- NA
        if(sunset == 'Inf') sunset <- NA
        
        # mean tidal derivative for metabolic day
        daytot <- mean(diff(x$Tide, na.rm = T))
        
        c(daytot, sunrise, sunset)
        
        }
      )
    names(tide_rngs) <- c('Date','daytot', 'sunrise', 'sunset')
    
    # get metab data from list
    dat_in <- met_ls[[x]]
    dat_in <- merge(dat_in, tide_rngs, by = 'Date', all.x = T)
    dat_in$month <- strftime(dat_in$Date, '%m')
    
    # Pg values correlated with tidal range during sunlight hours
    # Rt values correlated with tidal range during night hours
    # NEM values correlated with metabolic daily tidal range
    # NA if error is returned in correlation
    erf<- function(e) NA
    met_cor <- ddply(
        dat_in, 
        .variable = c('month'), 
        .fun = function(x){
          
          with(x, {c(
            Pg_obs = tryCatch(cor.test(Pg, sunrise)$estimate, error = erf),
            Rt_obs = tryCatch(cor.test(Rt, sunset)$estimate, error = erf),
            NEM_obs = tryCatch(cor.test(NEM, daytot)$estimate, error = erf),
            Pg_dtd = tryCatch(cor.test(Pg_dtd, sunrise)$estimate, error = erf),
            Rt_dtd =  tryCatch(cor.test(Rt_dtd, sunset)$estimate, error = erf),
            NEM_dtd = tryCatch(cor.test(NEM_dtd, daytot)$estimate, error = erf)
          )})
          
        }
      )
    names(met_cor) <- gsub('\\.cor$', '', names(met_cor))
    
    # DO and metab corrs combined
    res_sum <- data.frame(met_cor, do_obs = do_obs$cor, do_dtd = do_dtd$cor)
    
    # add column for sim name
    res_sum$L1 <- nm
    
    res_sum
    
  })
names(cor_res) <- case_regs
cor_res <- do.call('rbind', cor_res)
row.names(cor_res) <- 1:nrow(cor_res)

cor_res$site <- gsub('_wtreg_[0-9]*$', '', cor_res$L1)
cor_res$wins <- as.numeric(gsub('^.*_wtreg_', '', cor_res$L1))

# merge with case_grds
case_grds$wins <- as.numeric(row.names(case_grds))
cor_res <- merge(cor_res, case_grds, by = 'wins', all.x = T)

cor_res <- melt(cor_res, measure.var = c('Pg_obs', 'Rt_obs', 'NEM_obs', 'Pg_dtd', 'Rt_dtd', 'NEM_dtd', 'do_obs', 'do_dtd'))

# create columns for variable (DO, flux, etc.) and sub variable (obs, dtd)
cor_res$sub_var <- gsub('^.*_', '', cor_res$variable)
cor_res$var <- gsub('_.*$', '', cor_res$var)

# subset cor_res for window of 30 days, tidal prop 1
load('case_grds.RData')

pdf('cor_wins.pdf', height = 5.5, width = 9, family = 'serif')

for(row in 1:nrow(case_grds)){
  
  # plot
  to_plo <- cor_res[cor_res$wins == row, ]
  
  # reassign factor labels
  to_plo$month <- factor(to_plo$month, levels = c('01', '02', '03', '04',
    '05', '06', '07', '08', '09', '10', '11', '12'), 
    labels = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11',' 12')
    )
  to_plo$var <- factor(to_plo$var, levels = c('Pg', 'Rt', 'NEM', 'do'), 
    labels = c('Pg', 'Rt', 'NEM', 'DO'))
  to_plo$sub_var <- factor(to_plo$sub_var, levels = c('dtd', 'obs'), 
    labels = c('Filtered', 'Observed'))
  to_plo$site <- factor(to_plo$site, 
    levels = c('ELKVM', 'PDBBY', 'RKBMB', 'SAPDC'),
    labels = c('Elkhorn Slough', 'Padilla Bay', 'Rookery Bay', 
      'Sapelo Island')
    )
  
  tit_val <- case_grds[row, ]
  tit_val <- melt(tit_val)
  tit_val <- paste(paste(tit_val[, 1], tit_val[, 2]), collapse = ', ')
  p <- ggplot(to_plo, aes(x = factor(month), y = value, group = sub_var, colour = sub_var)) + 
    geom_line() +
    geom_point() + 
    geom_hline(yintercept = 0, linetype = 'dashed') + 
    facet_grid(site ~ var) + 
    scale_y_continuous(limits = c(-1, 1)) +
    ggtitle(tit_val)
    ylab('Correlation with tide') 
  
  print(p)
}
dev.off()

##
# heat map of correlations
# reassign factor labels for windows

to_plo <- cor_res

# average correlations by month
to_plo <- ddply(to_plo,
  .variable = c('wins', 'sub_var', 'var', 'site', 'dec_time', 'hour', 'Tide', 'L1', 'variable'), 
  .fun = function(x) mean(x$value, na.rm = T)
)

# remove observed data
to_plo <- to_plo[to_plo$sub_var == 'dtd', ]
to_plo$sub_var <- NULL

labs <- paste('Tide', unique(to_plo$Tide))
to_plo$Tide <- factor(to_plo$Tide, labels = labs)

mat_theme <- theme(
    panel.margin = unit(0, 'lines'), 
    strip.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.position = 'top',
    axis.text = element_text(size = 8)
    )

p1 <- ggplot(to_plo[to_plo$site == 'ELKVM', ], aes(x = factor(dec_time), y = factor(hour), z = V1, fill = V1)) + 
  geom_tile() +
  facet_grid(var ~ Tide) +
  scale_fill_gradientn(name = 'Correlation', 
    colours = brewer.pal(11, 'Spectral')#, limits = c(0.3, 0.7),
#     breaks = seq(0.3, 0.7, by = 0.1)
    ) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  xlab('Days') +
  ylab('Hours') + 
  theme_bw() +
  mat_theme


grid.arrange(p1, p2, ncol = 2)