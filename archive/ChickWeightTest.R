
#  from:https://medium.com/@marc.jacobs012/mixed-models-of-chicks-weight-a41413c99b51 




rm(list = ls())
# load libraries
library(haven)
library(lme4)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(sjmisc)
library(ggplot2)
theme_set(theme_sjplot())
library(arm)
library(ggplot2)
library(pbkrtest)
library(lmerTest)
library(plyr)
library(ggplot2)
library(kml)
library(effects)
library(lme4)
library(AICcmodavg)
library(parallel)
library(broom.mixed)
library(rstanarm)
library(tidybayes)
library(modelr)
library(RColorBrewer)
library(bayesplot)
library(MCMCglmm)
library(lme4)
library(ggstatsplot)
library(brms)
set.seed(123)

mod1 <- lm(weight ~ Time*Diet, 
           data = ChickWeight)
mod2 <- lmer(weight ~ Time*Diet +(1|Chick), 
             data = ChickWeight)
combine_plots(
  plotlist = list(
    ggcoefstats(mod1) +
      ggplot2::labs(x = parse(text = "'Regression coefficient' ~italic(beta)")),
    ggcoefstats(mod2) +
      ggplot2::labs(
        x = parse(text = "'Regression coefficient' ~italic(beta)"),
        y = "fixed effects"
      )),
  plotgrid.args = list(nrow = 2),
  annotation.args = list(title = "Relationship time + diet on weight of chickens"))


#### CHICKWEIGHT ####
dim(ChickWeight)
skimr::skim(ChickWeight)

ggplot(ChickWeight, 
       aes(x=as.factor(Time), y=weight, 
           fill=as.factor(Diet))) +
  geom_boxplot() +
  theme_bw() +
  ggtitle(paste("Boxplot of BW by time")) +
  scale_y_continuous("BW") +
  scale_x_discrete("Time") +
  labs(fill="Diet")+
  theme(text=element_text(size=10), 
        axis.title.x = element_text(), panel.grid.minor = element_blank())

ggplot(ChickWeight, aes(x=as.factor(Time), 
                        y=weight, 
                        group=Chick)) + 
  geom_point(alpha=0.6) + 
  geom_smooth(aes(group=Diet, colour=as.factor(Diet)), method="loess", size=1.5, se=F, alpha=0.8, span=0.5) +
  scale_y_continuous("BW") +
  scale_x_discrete("Time")+
  ggtitle(paste("Smoothed curve BW"))+
  labs(color="Diet")+
  theme_bw()+
  theme(text=element_text(size=10), axis.title.x = element_text(), 
        panel.grid.minor = element_blank())

ggplot(ChickWeight, aes(x=as.factor(Time), y=weight, group=Chick, colour=as.factor(Diet))) + 
  geom_point(alpha=0.6) + 
  geom_smooth(aes(group=1), size=1.5, se=F, alpha=0.8, span=0.5) +
  facet_grid(~Diet)+
  theme_bw()

summary(regr<-lm(weight~Time*Diet, data=ChickWeight))
par(mfrow=c(2,2))
plot(regr)
plot_model(regr)
plot(allEffects(regr))
dev.off()
model.i<-lmer(weight~1 + (1|Chick), ChickWeight, REML=FALSE) # random intercept with fixed mean - used to estimate intraclass correlation
summary(model.i) # fixed effect is the GRAND MEAN, which is 10.24 kg BW
plot(ChickWeight$Time, ChickWeight$weight, type="l"); abline(h=121.007, col="red", lwd=4)# shows how the grand mean does not in the very least capture the data
performance::icc(model.i)


fit<-lm(weight~1 + Time + factor(Chick, ordered=FALSE),
        ChickWeight)
summary(fit)
par(mfrow = c(2, 2))
plot(fit)
dev.off()
model.l<-lmer(weight~1 + Time + (1|Chick), ChickWeight, REML=FALSE) # random intercept + fixed time
summary(model.l) 
regr<-lm(ChickWeight$weight~ChickWeight$Time) 
plot(ChickWeight$Time, ChickWeight$weight, type="l"); abline(regr, col="red", lwd=4) # captures the effect of time better then the grand mean does

fit_augmented <- augment(model.l)
ggplot(fit_augmented, aes(x=Time, group=as.factor(Chick)))+
  geom_line(aes(y=weight), colour="black")+
  geom_line(aes(y=.fitted), colour="red", alpha=1)+
  xlab("Time")+ylab("BW")+
  ggtitle("Observed vs Predicted")+
  theme_bw()+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(fit_augmented, aes(x=weight, y=fitted(model.l), colour=weight-fitted(model.l)))+
  geom_point()+scale_colour_gradient2(low="red",mid="blue",high="red")+
  geom_abline(intercept=0, slope=1, col="black", linetype="dashed")+
  labs(title = "Calibration", y="Predicted", x="Observed", colour="Difference")+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

model.q<-lmer(weight~1 + Time+ I(Time^2) + (1|Chick), ChickWeight, REML=FALSE)
summary(model.q) # BW=7.6313063+0.1607967*time+0.0056027*time^2
regr2<-lm(weight~Time+I(Time^2), data=ChickWeight) 
plot(ChickWeight$Time, ChickWeight$weight, type="l");
timevalues<-seq(0,21,1); 
predicted<-predict(regr2, list(Time=timevalues, 
                               Time2=I(timevalues^2))); 
lines(timevalues, predicted, col="red", lwd=4) 
abline(regr, col="blue", lwd=4) 

mynames<-c("I", "L", "Q")
myaicc<-as.data.frame(aictab(cand.set=list(model.i,model.l,model.q), modnames=mynames))[, -c(5,7)]
myaicc$eratio<-max(myaicc$AICcWt)/myaicc$AICcWt
data.frame(Model=myaicc[,1],round(myaicc[, 2:7],4)) # prefers quadratic, but all it rules out is the intercept model


model.l.t.r<-lmer(weight~Time + Diet + (Time|Chick), ChickWeight, REML=FALSE) # random intercept & slope time + effect of treatment
summary(model.l.t.r)

fit_augmented <- augment(model.l.t.r)
ggplot(fit_augmented, aes(x=Time, 
                          group=as.factor(Chick)))+
  geom_line(aes(y=weight), colour="black")+
  geom_line(aes(y=.fitted), colour="red", alpha=1)+
  xlab("Time")+ylab("BW")+
  ggtitle("Observed vs Predicted")+
  theme_bw()+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(fit_augmented, aes(x=weight, 
                          y=fitted(model.l.t.r), 
                          colour=weight-fitted(model.l.t.r)))+
  geom_point()+scale_colour_gradient2(low="red",mid="blue",high="red")+
  geom_abline(intercept=0, slope=1, col="black", linetype="dashed")+
  labs(title = "Calibration", y="Predicted", x="Observed", colour="Difference")+
  theme_bw()+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Extracting fixed and random effects
head(ranef(model.l.t.r))
head(fixef(model.l.t.r))
# Extracting covariance matrix
round(vcov(model.l.t.r),3)
plot_model(model.l.t.r, grid=T, sort.est="sort.all", y.offset=1.0)
plot_model(model.l.t.r, grid=T, type="diag")
plot_model(model.l.t.r, type="eff")
plot_model(model.l.t.r, show.values = TRUE, value.offset = .3)
plot_model(model.l.t.r, type="re") 
plot_model(model.l.t.r, type="std")  
plot_model(model.l.t.r, type="slope") 
plot_model(model.l.t.r, type="resid")
plot_model(model.l.t.r, type="pred", grid=T, terms=c("Time [all]", "Diet"))
plot_model(model.l.t.r, type="pred", pred.type = "re", grid=T, terms=c("Time [all]", "Diet"))
plot_model(model.l.t.r, type="pred", show.data = TRUE, pred.type = "fe", grid=T, 
           terms=c("Time [all]", "Diet"))
tab_model(model.l.t.r, 
          show.std=T, 
          show.r2=T, 
          show.aic=T, 
          show.aicc = T, 
          show.dev=T, 
          show.re.var=T)

model.l.t.r<-lmer(weight~Time + Diet + (Time|Chick), ChickWeight, REML=FALSE)
summary(model.l.t.r)
model.l.t.r.bayes<-stan_lmer(weight~Time + Diet + (Time|Chick), 
                             data=ChickWeight)
summary(model.l.t.r.bayes)
plot(model.l.t.r.bayes)
prior_summary(object = model.l.t.r.bayes)
summary(model.l.t.r.bayes, 
        probs = c(0.025, 0.975),
        digits = 2)
plot(model.l.t.r.bayes, "rhat")
plot(model.l.t.r.bayes, "ess")

sims <- as.matrix(model.l.t.r.bayes)
para_name <- colnames(sims)
mu_a_sims <- as.matrix(model.l.t.r.bayes, 
                       pars = "(Intercept)")
u_sims <- as.matrix(model.l.t.r.bayes, 
                    regex_pars = "b\\[\\(Intercept\\) Chick\\:")
a_sims <- as.numeric(mu_a_sims) + u_sims          
s_y_sims <- as.matrix(model.l.t.r.bayes, 
                      pars = "sigma")
s__alpha_sims <- as.matrix(model.l.t.r.bayes, 
                           pars = "Sigma[Chick:(Intercept),(Intercept)]")
a_mean <- apply(X = a_sims,     # posterior mean
                MARGIN = 2,
                FUN = mean)
a_sd <- apply(X = a_sims,       # posterior SD
              MARGIN = 2,
              FUN = sd)
a_quant <- apply(X = a_sims, 
                 MARGIN = 2, 
                 FUN = quantile, 
                 probs = c(0.025, 0.50, 0.975))
a_quant <- data.frame(t(a_quant))
names(a_quant) <- c("Q2.5", "Q50", "Q97.5")
a_df <- data.frame(a_mean, a_sd, a_quant)
round(head(a_df), 2)
a_df <- a_df[order(a_df$a_mean), ]
a_df$a_rank <- c(1 : dim(a_df)[1])  
ggplot(data = a_df, 
       aes(x = a_rank, 
           y = a_mean)) +
  geom_pointrange(aes(ymin = Q2.5, 
                      ymax = Q97.5),
                  position = position_jitter(width = 0.1, 
                                             height = 0)) + 
  geom_hline(yintercept = mean(a_df$a_mean), 
             size = 0.5, 
             col = "red") + 
  scale_x_continuous("Rank", 
                     breaks = seq(from = 0, 
                                  to = 80, 
                                  by = 5)) + 
  scale_y_continuous(expression(paste("varying intercept, ", alpha[j]))) + 
  theme_bw( base_family = "serif")

sims <- as.matrix(model.l.t.r.bayes)
para_name <- colnames(sims)
mu_a_sims <- as.matrix(model.l.t.r.bayes, 
                       pars = "Time")
u_sims <- as.matrix(model.l.t.r.bayes, 
                    regex_pars = "b\\[\\Time\\ Chick\\:")
a_sims <- as.numeric(mu_a_sims) + u_sims          
s_y_sims <- as.matrix(model.l.t.r.bayes, 
                      pars = "sigma")
s__alpha_sims <- as.matrix(model.l.t.r.bayes, 
                           pars = "Sigma[Chick:Time,Time]")
a_mean <- apply(X = a_sims,     # posterior mean
                MARGIN = 2,
                FUN = mean)
a_sd <- apply(X = a_sims,       # posterior SD
              MARGIN = 2,
              FUN = sd)
a_quant <- apply(X = a_sims, 
                 MARGIN = 2, 
                 FUN = quantile, 
                 probs = c(0.025, 0.50, 0.975))
a_quant <- data.frame(t(a_quant))
names(a_quant) <- c("Q2.5", "Q50", "Q97.5")
a_df <- data.frame(a_mean, a_sd, a_quant)
round(head(a_df), 2)
a_df <- a_df[order(a_df$a_mean), ]
a_df$a_rank <- c(1 : dim(a_df)[1])  
ggplot(data = a_df, 
       aes(x = a_rank, 
           y = a_mean)) +
  geom_pointrange(aes(ymin = Q2.5, 
                      ymax = Q97.5),
                  position = position_jitter(width = 0.1, 
                                             height = 0)) + 
  geom_hline(yintercept = mean(a_df$a_mean), 
             size = 0.5, 
             col = "red") + 
  scale_x_continuous("Rank", 
                     breaks = seq(from = 0, 
                                  to = 80, 
                                  by = 5)) + 
  scale_y_continuous(expression(paste("varying slope, ", alpha[j]))) + 
  theme_bw( base_family = "serif")

ChickWeight %>%
  add_epred_draws(model.l.t.r.bayes, 
                  ndraws=50) %>%
  ggplot(aes(x = Time, 
             y = weight, 
             color = as.factor(Diet), 
             group=Chick)) +
  geom_line(aes(y = .epred, group = paste(Chick, .draw)), alpha = 0.25) +
  geom_point(data = ChickWeight, color="black")+
  facet_wrap(~Diet)+
  labs(color="Diet")+
  theme_bw()

model.l.t.r.bayes2<-brms::brm(weight~Time + Diet + (Time|Chick), 
                              data=ChickWeight)
summary(model.l.t.r.bayes2)
plot(model.l.t.r.bayes2)

posterior <- as.matrix(model.l.t.r.bayes2)
mcmc_areas(posterior,
           pars = c("b_Intercept", "b_Time", "b_Diet2","b_Diet3","b_Diet4"),
           prob = 0.8) 

pp_check(model.l.t.r.bayes2, ndraws = 500)
pp_check(model.l.t.r.bayes2, type = "error_hist", ndraws = 11)
pp_check(model.l.t.r.bayes2, type = "scatter_avg", ndraws = 100)
pp_check(model.l.t.r.bayes2, type = "stat_2d")
pp_check(model.l.t.r.bayes2, type = "loo_pit")

model.l.t.r.bayes2 %>%
  posterior_predict(draws = 500) %>%
  ppc_stat_grouped(y = ChickWeight$weight,
                   group = ChickWeight$Diet,
                   stat = "median")

np <- nuts_params(model.l.t.r.bayes2 )
mcmc_nuts_energy(np) + 
  ggtitle("NUTS Energy Diagnostic")

brms::get_prior(weight~Time + Diet + (Time|Chick), 
                data=ChickWeight)
model.l.t.r.bayes2<-brms::brm(weight~Time + 
                                Diet + 
                                (Time|Chick), 
                              data=ChickWeight,
                              prior = c(
                                prior(normal(100, 20),   class = Intercept),
                                prior(normal(20, 10),   class = b,  coef=Diet2),
                                prior(normal(10, 5),   class = b,  coef=Diet3),
                                prior(normal(40, 10),   class = b,  coef=Diet4),
                                prior(gamma(2, 2), class = sd, coef=Intercept, group=Chick),
                                prior(gamma(2, 2), class = sd, coef=Time,      group=Chick),
                                prior(gamma(2,2), class = sigma)),
                              chains=4, cores=6)
summary(model.l.t.r.bayes2)
plot(model.l.t.r.bayes2)

pp_check(model.l.t.r.bayes2, ndraws = 500)
pp_check(model.l.t.r.bayes2, type = "error_hist", ndraws = 11)
pp_check(model.l.t.r.bayes2, type = "scatter_avg", ndraws = 100)
pp_check(model.l.t.r.bayes2, type = "stat_2d")
pp_check(model.l.t.r.bayes2, type = "loo_pit")

ppc_intervals(
  y = ChickWeight$weight,
  yrep = posterior_predict(model.l.t.r.bayes2),
  x = ChickWeight$weight,
  prob = 0.5) +
  labs(
    x = "Time",
    y = "Weight",
    title = "50% posterior predictive intervals") +
  panel_bg(fill = "gray95", color = NA) +
  grid_lines(color = "white")+theme_bw()

ppc_intervals_grouped(
  y = ChickWeight$weight,
  yrep = posterior_predict(model.l.t.r.bayes2),
  x = ChickWeight$weight,
  group=ChickWeight$Diet,
  prob = 0.5) +
  labs(
    x = "Time",
    y = "Weight",
    title = "50% posterior predictive intervals") +
  panel_bg(fill = "gray95", color = NA) +
  grid_lines(color = "white")+theme_bw()

ppc_ribbon_grouped(
  y = ChickWeight$weight,
  yrep = posterior_predict(model.l.t.r.bayes2),
  x = ChickWeight$weight,
  group=ChickWeight$Diet,
  prob = 0.5) +
  labs(
    x = "Time",
    y = "Weight",
    title = "50% posterior predictive intervals") +
  panel_bg(fill = "gray95", color = NA) +
  grid_lines(color = "white")+theme_bw()

ChickWeight %>%
  add_epred_draws(model.l.t.r.bayes2, 
                  ndraws=50) %>%
  ggplot(aes(x = Time, 
             y = weight, 
             color = as.factor(Diet), 
             group=Chick)) +
  geom_line(aes(y = .epred, group = paste(Chick, .draw)), alpha = 0.25) +
  geom_point(data = ChickWeight, color="black")+
  facet_wrap(~Diet)+
  labs(color="Diet")+
  theme_bw()


ChickWeight %>%
  dplyr::group_by(Diet) %>%
  add_epred_draws(model.l.t.r.bayes2, ndraws = 20) %>%
  ggplot(aes(x = Time, y = weight, color = ordered(Diet))) +
  geom_line(aes(y = .epred, group = paste(Diet,Chick, .draw))) +
  geom_point(data = ChickWeight) +
  scale_color_brewer(palette = "Dark2") + theme_bw()

ChickWeight %>%
  dplyr::group_by(Diet) %>%
  add_predicted_draws(model.l.t.r.bayes2, ndraws = 20) %>%
  ggplot(aes(x = Time, y = weight, color = ordered(Diet),fill = ordered(Diet))) +
  stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50), alpha = 1/4) +
  geom_point(data = ChickWeight) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Dark2")+theme_bw()

ChickWeight %>%
  dplyr::group_by(Diet) %>%
  add_predicted_draws(model.l.t.r.bayes2) %>%
  ggplot(aes(x = Time, y = weight)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.99, .95, .8, .5), 
                  color = brewer.pal(5, "Blues")[[5]]) +
  geom_point(data = ChickWeight) +
  scale_fill_brewer() +
  facet_grid(. ~ Diet, space = "free_x", scales = "free_x")+theme_bw()

grid = ChickWeight %>%data_grid(Diet, Time, Chick)
means = grid %>%add_epred_draws(model.l.t.r.bayes2)
preds = grid %>%add_predicted_draws(model.l.t.r.bayes2)
ChickWeight %>%
  ggplot(aes(x = weight, y = Diet)) +
  stat_halfeye(aes(x = .epred), scale = 0.6, 
               position = position_nudge(y = 0.175), data = means) +
  stat_interval(aes(x = .prediction), data = preds) +
  geom_point(data = ChickWeight) +
  scale_color_brewer()+theme_bw()

model.l.t.r.bayes2%>%
  recover_types(ChickWeight) %>%
  emmeans::emmeans( ~ Diet, data = ChickWeight) %>%
  emmeans::contrast(method = "pairwise") %>%
  gather_emmeans_draws() %>%
  ggplot(aes(x = .value, y = contrast)) +
  stat_halfeye()+theme_bw()


model.l.t.r.bayes2%>%
  recover_types(ChickWeight) %>%
  emmeans::emmeans( ~ Diet | Time, data = ChickWeight) %>%
  emmeans::contrast(method = "pairwise") %>%
  gather_emmeans_draws() %>%
  ggplot(aes(x = .value, y = contrast)) +
  stat_halfeye()+theme_bw()

get_prior(weight~Time*Diet + (Time|Chick), 
          data=ChickWeight)
model.l.t.r.bayes3<-brms::brm(weight~Time*Diet + (Time|Chick),
                              prior = c(
                                prior(normal(100, 20), class = Intercept),
                                prior(normal(3, 0.5),  class = b,  coef=Time),
                                prior(normal(20, 10),  class = b,  coef=Diet2),
                                prior(normal(10, 5),   class = b,  coef=Diet3),
                                prior(normal(40, 10),  class = b,  coef=Diet4),
                                prior(normal(3, 1),  class = b,  coef=Time:Diet2),
                                prior(normal(-2, 1),   class = b,  coef=Time:Diet3),
                                prior(normal(8, 2),  class = b,  coef=Time:Diet4),
                                prior(gamma(2, 2), class = sd, coef=Intercept, group=Chick),
                                prior(gamma(2, 2), class = sd, coef=Time,      group=Chick),
                                prior(gamma(2,2), class = sigma)),
                              chains=4, cores=6,
                              data=ChickWeight)
model.l.t.r.bayes3$prior
plot(model.l.t.r.bayes3)
summary(model.l.t.r.bayes3)

pp_check(model.l.t.r.bayes3, ndraws = 500)

ppc_intervals_grouped(
  y = ChickWeight$weight,
  yrep = posterior_predict(model.l.t.r.bayes3),
  x = ChickWeight$weight,
  group=ChickWeight$Diet,
  prob = 0.5) +
  labs(
    x = "Time",
    y = "Weight",
    title = "50% posterior predictive intervals") +
  panel_bg(fill = "gray95", color = NA) +
  grid_lines(color = "white")+theme_bw()

ppc_intervals_grouped(
  y = ChickWeight$weight,
  yrep = posterior_predict(model.l.t.r.bayes3),
  x = ChickWeight$weight,
  group=ChickWeight$Time,
  prob = 0.5) +
  labs(
    x = "Time",
    y = "Weight",
    title = "50% posterior predictive intervals") +
  panel_bg(fill = "gray95", color = NA) +
  grid_lines(color = "white")+theme_bw()


pp_check(model.l.t.r.bayes3, type = "error_hist", ndraws = 11)
pp_check(model.l.t.r.bayes3, type = "scatter_avg", ndraws = 100)
pp_check(model.l.t.r.bayes3, type = "stat_2d")
pp_check(model.l.t.r.bayes3, type = "loo_pit")


model.l.t.r.bayes3%>%
  recover_types(ChickWeight) %>%
  emmeans::emmeans( ~ Diet | Time, data = ChickWeight) %>%
  emmeans::contrast(method = "pairwise") %>%
  gather_emmeans_draws() %>%
  ggplot(aes(x = .value, y = contrast)) +
  facet_wrap(~Time)+
  stat_halfeye()+theme_bw()




view rawChickweight Bayesian hosted with ❤ by GitHub