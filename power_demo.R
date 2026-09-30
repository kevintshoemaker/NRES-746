
rm(list=ls())

## Set global parameters --------------

base_params = list()

base_params$N_0  = 1000
base_params$p = 0.02
base_params$k = 3
base_params$s = 1
base_params$y = 25
mindecline_25 = 0.75
base_params$alpha = 0.05
base_params$lam = mindecline_25^(1/base_params$y)


## custom functions ----------


p_bout <- function(p,k)  1 - (1-p)^k

   # tests
p_bout(0.02,3)
1- 0.98^3

true_pop_traj = function(n0, l, y) n0 *l^(1:y)


true_pop_traj(100,1,10)

plot(1:50,true_pop_traj(100,0.95,50))

survyears <- function(y,s) seq(1,y,by=s)
survyears(10,2)

pars=base_params
sim_counts = function(pars){
  ys = survyears(pars$y,pars$s)
  traj = true_pop_traj(pars$N_0,pars$lam,pars$y)
  rbinom(length(ys),round(traj[ys]),p_bout(pars$p,pars$k))
}

p_test = base_params; p_test$lam = 0.95
sim_counts(p_test)


counts=sim_counts(pars); years=survyears(pars$y,pars$s)
detect_decline = function(years,counts,a){
  mod = lm(log(counts)~years)
  is_decline = unname(coef(mod)[2] < 0)
  is_sig = summary(mod)$coefficients[2,4] < a
  is_decline && is_sig 
}

counts=sim_counts(pars); years=survyears(pars$y,pars$s)
detect_decline(counts,years,0.05)

do_rep = function(pars){
  counts=sim_counts(pars); years=survyears(pars$y,pars$s)
  detect_decline(counts,years,pars$alpha)
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
  sy = survyears(params$y,params$s)
  grid$cost[g] = length(sy)*(2000+params$k*200)
  grid$pow[g] = get_power(1000,params)
}

grid2 = grid[grid$pow > 0.75,]

grid2 = grid2[order(grid2$cost),]





