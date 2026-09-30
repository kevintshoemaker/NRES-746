
## Set global parameters --------------

base_params = list()

base_params$N_0  = 1000
base_params$p = 0.02
base_params$k = 3
base_params$s = 1
base_params$y = 25
base_params$mindecline_25 = 0.75
base_params$alpha = 0.05
base_params$lam = base_params$mindecline_25^(1/base_params$y)


## custom functions ----------


p_bout <- function(pars){
  1 - (1-pars$p)^pars$k
}

   # tests
p_test = base_params; p_test$k = Inf   # try known edge cases
p_bout(p_test)
1- 0.98^3


true_pop_traj = function(pars){
  pars$N_0 * pars$lam^(1:pars$y)
}

p_test = base_params; p_test$lam = 0.98
true_pop_traj(p_test)


survyears <- function(pars) seq(1,pars$y,by=pars$s)

sim_counts = function(pars){
  rbinom(length(survyears(pars)),round(true_pop_traj(pars)[survyears(pars)]),p_bout(pars))
}

sim_counts(p_test)

counts=sim_counts(pars); years=survyears(pars)
detect_decline = function(years,counts,pars){
  mod = lm(log(counts)~years)
  is_decline = unname(coef(mod)[2] < 0)
  is_sig = summary(mod)$coefficients[2,4] < pars$alpha
  is_decline && is_sig 
}

counts=sim_counts(pars); years=survyears(pars)
detect_decline(counts,years,p_test)

do_rep = function(pars){
  counts=sim_counts(pars); years=survyears(pars)
  detect_decline(counts,years,pars)
}

get_power = function(reps,pars){
  sum(replicate(reps,do_rep(pars)))/reps
}

get_power(1000,pars)

### set up test grid

grid = expand.grid(k=1:15, s=1:5)
grid$cost = NA
grid$pow = NA

g=1
for(g in 1:nrow(grid)){
  params = base_params
  params$k = grid$k[g]
  params$s = grid$s[g]
  sy = survyears(params)
  grid$cost[g] = length(sy)*(2000+params$k*200)
  grid$pow[g] = get_power(1000,params)
}

grid2 = grid[grid$pow > 0.75,]

grid2 = grid2[order(grid2$cost),]





