
#  NRES 746, Lecture 2  ---------------------------           
##  University of Nevada, Reno     
##  Working with Probabilities                   


# Classic Urn Example ------------------------

Urn_summary <- c(
  red=104,
  blue=55,
  green=30
)

set.seed(1001)
Urn_contents <- sample(rep(names(Urn_summary),times=Urn_summary))

table(Urn_contents)   # verify that the contents are what you expect!


P_blue <- sum(Urn_contents=="blue")/length(Urn_contents) # probability of drawing a blue sphere
P_blue


Prob <- table(Urn_contents)/length(Urn_contents)    # probability of drawing each type of sphere
Prob


as.numeric( Prob["blue"] + Prob["red"] )     # probability of drawing a blue or red sphere


as.numeric( Prob["blue"] + Prob["red"] + Prob["green"] )      # P(blue OR green)


# Question: What is the probability of drawing a blue **AND THEN** a red sphere? 

#[your command here]    # P(blue AND THEN red)


#### Question: What is the probability of drawing a blue and a red sphere in two consecutive draws (but in no particular order)? 

#[your command here]    # P(blue AND THEN red)


# Urn example #2: color and shape --------------------------

Urn_summary <- matrix(NA,2,2,dimnames = list(Color=c("red","blue"),Shape=c("sphere","cube") ))
Urn_summary["red","sphere"] <- 39   # contents of new urn
Urn_summary["blue","sphere"] <- 76
Urn_summary["red","cube"] <- 101   
Urn_summary["blue","cube"] <- 25
Urn_summary

Prob = Urn_summary/sum(Urn_summary)


Prob_Shape <- colSums(Urn_summary)/sum(Urn_summary)  # marginal probabilities of shape
Prob_Shape

Prob_Color <- rowSums(Urn_summary)/sum(Urn_summary)    # marginal probabilities of color
Prob_Color


as.numeric( Prob_Color["blue"] * Prob_Shape["cube"])      # joint probability of drawing a blue object that is a cube

## NOTE: if the above answer is not correct, please correct it! 


as.numeric( Prob_Color["blue"] + Prob_Shape["cube"])        # probability of drawing something blue or something cube-shaped...

## NOTE: if the above answer is not correct, please correct it!  



Prob["blue","cube"] / Prob_Shape["cube"]   # probability of drawing a blue object, given it is a cube


as.numeric( (Prob["blue","cube"] / Prob_Shape["cube"]) * Prob_Shape["cube"])   # probability of drawing a blue cube... using conditional probabilities

Prob["blue","cube"]   # check answer to make sure it's right


# unconditional probability of drawing a blue item.  Seems too complicated, but this method of computing unconditional probabilities will prove useful as we get into Bayesian statistics!
uncond_prob_blue <- (Prob["blue","cube"] /  Prob_Shape["cube"]) * Prob_Shape["cube"] + 
              (Prob["blue","sphere"] / Prob_Shape["sphere"]) * Prob_Shape["sphere"]       

as.numeric(uncond_prob_blue)


Prob_Shape <- colSums(Urn_summary)/sum(Urn_summary)  # marginal probabilities of shape
Prob_Shape

Prob_Color <- rowSums(Urn_summary)/sum(Urn_summary)    # marginal probabilities of color
Prob_Color

Prob_Color["blue"]      # marginal probability of drawing a blue object (across all possible shapes)


# Medical example (positive predictive value) --------------------------

Prob_Disease <- c(1e-4, 1-1e-4)     # marginal probability of disease vs no disease
names(Prob_Disease) <- c("yes","no")                # make it a named vector!
Prob_Disease


## compute the unconditional probability of testing positive

as.numeric( 1*Prob_Disease["yes"] + 0.01*Prob_Disease["no"] )    # Prob(+test|Disease)*Prob(Disease) + Prob(+test|no Disease)*Prob(no Disease)


# Monty Hall simulation code --------------
#   (code by Corey Chivers 2012)
 
monty<-function(strat='stay',N=1000,print_games=TRUE){
  doors<-1:3 #initialize the doors behind one of which is a good prize
  win<-0 #to keep track of number of wins
  
  for(i in 1:N){
    prize<-floor(runif(1,1,4)) #randomize which door has the good prize
    guess<-floor(runif(1,1,4)) #guess a door at random
    
    ## Reveal one of the doors you didn't pick which has a bum prize
    if(prize!=guess)
      reveal<-doors[-c(prize,guess)]
    else
      reveal<-sample(doors[-c(prize,guess)],1)
    
    ## Stay with your initial guess or switch
    if(strat=='switch')
      select<-doors[-c(reveal,guess)]
    if(strat=='stay')
      select<-guess
    if(strat=='random')
      select<-sample(doors[-reveal],1)
    
    ## Count up your wins
    if(select==prize){
      win<-win+1
      outcome<-'Winner!'
    }else
      outcome<-'Loser!'
    
    if(print_games)
      cat(paste('Guess: ',guess,
          '\nRevealed: ',reveal,
          '\nSelection: ',select,
          '\nPrize door: ',prize,
          '\n',outcome,'\n\n',sep=''))
  }
  cat(paste('Using the ',strat,' strategy, your win percentage was ',win/N*100,'%\n',sep='')) #Print the win percentage of your strategy
}


# run the monty hall code!

monty(strat="stay",print_games=FALSE)


# run the monty hall code!

monty(strat="switch",print_games=FALSE)


# Probability distributions in R  ---------------------

mean <- 5
rpois(10,mean)    # the random numbers are integers with no decimal component

## Discrete -------------------------

             # plot a discrete distribution!
xvals <- seq(0,15,1)
probs <- dpois(xvals,lambda=mean)
names(probs) <- xvals
               
barplot(probs,ylab="Probability",main="Poisson distribution (discrete)")

barplot(cumsum(probs),ylab="Cumulative Probability",main="Poisson distribution (discrete)")   # cumulative distribution

sum(probs)   # just to make sure it sums to 1!  Does it??? 


## Continuous  --------------------

alpha = 0.5
beta = 0.5

rbeta(10,alpha,beta)

curve(dbeta(x,alpha,beta))   # probability density

curve(pbeta(x,alpha,beta))   # cumulative distribution

integrate(f=dbeta,lower=0,upper=1,shape1=alpha,shape2=beta)    # just to make sure it integrates to 1!!


## Binomial -----------------

size <- 10
prob <- 0.3
rbinom(10,size,prob)

xvals <- seq(0,size,1)
probs <- dbinom(xvals,size,prob)
names(probs) <- xvals
               
barplot(probs,ylab="Probability",main="Binomial distribution")

barplot(cumsum(probs),ylab="Cumulative Probability",main="Binomial distribution")   # cumulative distribution

sum(probs)   # just to make sure it sums to 1!  Does it???


## Gaussian (Normal) ---------------

mean = 7.1
stdev = 1.9

rnorm(10,mean,stdev)

curve(dnorm(x,mean,stdev),0,15)   # probability density

curve(pnorm(x,mean,stdev),0,15)   # cumulative distribution

integrate(f=dnorm,lower=-Inf,upper=Inf,mean=mean,sd=stdev)    # just to make sure it integrates to 1!!

