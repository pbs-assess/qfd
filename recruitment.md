Recruitment exploraiton
================

##expnded Beverton-Holt (eq 3.6) where `a` and `b` are density
dependent, and independent additive effects on mortality, `Z`, and `T`
is age at recruitment.

``` r
S <- 1:100 #vector of spawners
T <- 4 #age @ recruitment
f <- 10  #avg net fecundity
a1 <- 1 #DI effect on mortality
b1 <- 0.1 #DD effect on mortality

R_L_BH <- S/((exp(a1*T)/f) + (b1/a1) + (exp(a1*T) - 1)*S)

plot(R_L_BH~S)
```

![](recruitment_files/figure-gfm/long%20beverton-holt-1.png)<!-- -->

We can get the more normal parameters, alpha and beta, with some math.
What do a1 and b1 *really mean*? How would one choose a reasonable
starting value for them for a population? Can you empirically estimate
these? Or they more abstract?

``` r
alpha <- f*exp(-a1*T)
beta <- f*b1*(1-exp(-a1*T))/a1
```

Then plug these values into the Beverton-Holt

``` r
R_BH <- (alpha*S)/(1+beta*S) #why TF is this different than Rec above?

plot(R_BH~S)
```

![](recruitment_files/figure-gfm/beverton%20holt-1.png)<!-- -->

Now we do the same as above by formulating the Ricker (eq. 3.8) with
additive DI and DD effects on mortality. **I’m going to try keeping the
a1 and b1 as above, but I’m not sure this is proper**  
- It’s probably not proper because in the book they use *a2* and *b2*,
meaning these could be the same thing but different numbers  
- and because `Z` us calculated a different way? **Need to revisit this
issue**

``` r
#S <- seq(0, 10, by = 0.1) # helper to toggle S around
R_L_R <- f*exp(-a1*T)*S*exp(-b1*T*S)

plot(R_L_R ~ S) #my code & math is right I just picked weird values for a and B I think... 
```

![](recruitment_files/figure-gfm/long%20Ricker-1.png)<!-- -->

If everything is right I *should be* able to write the Ricker with the
alpha and betas I calculated earlier?

``` r
R_R <- alpha*S*exp(-beta*S)

plot(R_R ~ S) #nope! - weird 
```

![](recruitment_files/figure-gfm/regular%20ricker-1.png)<!-- -->

Ludwig-Walters model, where the DD mortality term is a power function of
the Spawning stock `Zt = a2 + b2(S^gamma)`
