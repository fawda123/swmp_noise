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

load('met_ls_inst.RData')
load('met_ls.RData')

case <- 'PDBBY'

##
# all PAR
to.plo <- met_ls_inst[[grep(case, names(met_ls_inst))]]

p <- ggplot(to.plo, aes(DateTimeStamp, Tide, colour=TotPAR)) + 
  geom_line(size = 1) + 
  scale_x_datetime(name = '',expand=c(0,0)) + 
  scale_y_continuous('Depth (m)', expand=c(0,0), limits=c(0.5,4.5)) +
  scale_colour_gradient(low='blue', high='yellow',
    guide = guide_colorbar(direction = "horizontal",barheight= 0.3)) +
  theme_bw() +
  theme(legend.position = c(1,1), legend.justification = c(1,1),
    legend.background = element_rect(colour='black')) 
p

##
# subsets by case and date range
# plot metab and tide

dat.rng<-as.Date(c('2012-03-01','2012-03-27')) 

met_subs <- met_ls[[grep(case, names(met_ls_inst))]]
inst_subs <- met_ls_inst[[grep(case, names(met_ls_inst))]]

# subset data
met.rng <- met_subs$Date<=dat.rng[2] & met_subs$Date>=dat.rng[1]
met_subs <- met_subs[met.rng,]

inst_subs$Date <- as.Date(inst_subs$DateTimeStamp)
inst.rng <- inst_subs$Date<=dat.rng[2] & inst_subs$Date>=dat.rng[1]
inst_subs <- inst_subs[inst.rng,]

##
# custom theme, mod of theme_bw

my_theme <- theme(
  legend.title = element_blank(),legend.position = 'top',
  axis.title.x = element_blank(),legend.box= 'horizontal',
  plot.margin= unit(c(0, 1, 0, 1), "lines") # top right bottom left
  )

# function for setting range on y axis
rng.fun<-function(vec.in){
  rngs<-range(vec.in,na.rm=T)
  buffs<-0.07*abs(diff(rngs))
  c(rngs[1]-buffs,rngs[2]+buffs)
  }

##
# metab plot
to_plo1 <- melt(met_subs, id.var = c('Date'), 
  measure.var = grep('Pg|Rt|NEM', names(met_subs), value = T)
  )
to_plo1$Input <- 'Observed'
to_plo1$Input[grep('dtd', to_plo1$variable)] <- 'Detided'
to_plo1$Input <- factor(to_plo1$Input, levels = c('Observed', 'Detided'))
to_plo1$variable <- gsub('_dtd', '', to_plo1$variable)

ylab<-expression(paste('DO (g ',m^-2, d^-1, ')'))
p1 <- ggplot(to_plo1, aes(x = Date, y = 0.032 * value, group = variable,
    colour = variable)) +
  geom_line() +
  theme_bw() +
  geom_point(size = 2) +
  facet_wrap(~Input, ncol = 1) +
  scale_y_continuous(ylab)  +
  my_theme

##
# DO plot
to_plo2 <- met.day.fun(inst_subs, case)
names(to_plo2)[names(to_plo2) %in% 'variable'] <- 'solar'
ggpoly <- poly.fun(to_plo2$solar, to_plo2)

ylab<-expression(paste('DO (mg ',L^-1,')'))
p2 <- ggplot(to_plo2, aes(x = DateTimeStamp)) + 
  ggpoly +
  geom_line(aes(y = DO_obs, colour = 'Observed')) +
  geom_line(aes(y = DO_dtd, colour = 'Detided')) +
  coord_cartesian(ylim = rng.fun(to_plo2$DO_obs)) +
  scale_fill_manual(values='orange',labels='Day') +
  theme_bw() +
  scale_y_continuous(ylab)  +
  my_theme

## 
# tide plot
to_plo3 <- met.day.fun(inst_subs, case)
names(to_plo2)[names(to_plo2) %in% 'variable'] <- 'solar'
ggpoly <- poly.fun(to_plo2$solar, to_plo2)

ylab<-expression(paste('Height (m)'))
p3 <- ggplot(to_plo3, aes(x = DateTimeStamp)) + 
  ggpoly +
  geom_line(aes(y = Tide, colour = TotPAR), size = 1.2) +
  coord_cartesian(ylim = rng.fun(to_plo3$Tide)) +
  scale_fill_manual(values='orange',labels='Day') +
  theme_bw() +  
  my_theme + 
  scale_y_continuous(ylab) +
  scale_colour_gradient(low='blue', high='yellow',
    guide = guide_colorbar(direction = "horizontal",barheight= 0.3))

# Get the widths
pA <- ggplot_gtable(ggplot_build(p1))
pB <- ggplot_gtable(ggplot_build(p2))
pC <- ggplot_gtable(ggplot_build(p3))
maxWidth = unit.pmax(pA$widths[2:3], pB$widths[2:3], 
                     pC$widths[2:3])

# Set the widths
pA$widths[2:3] <- maxWidth
pB$widths[2:3] <- maxWidth
pC$widths[2:3] <- maxWidth

grid.arrange(pA, pB, pC, ncol = 1, heights = c(1.7, 1, 1))

