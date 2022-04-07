HowToInvadeTheOccupiedNiche\_Analyses
================
Mathias Wegner, AWI
07/01/2022

## R Markdown

This is an R Markdown of the analyses accompanying the manuscript
‘Invading the occupied niche: How a parasitic copepod of introduced
oysters can expel a congener from native mussels’ by Marieke E. Feis,
Leo Gottschalck, Lena Ruf, Franziska Theising, Felicitas Demann, K.
Mathias Wegner.

Functions:

``` r
#defining standard error
se = function(x)
{
    sd(x, na.rm=T)/sqrt(length(x))
}

#perform phase shift analysis
# parameters:
# x, y: numeric vectors for shifted correlation calucation
# left_shift, right_shift: integers for how many intervals series will be pushed left and right
# value:
# a list 1) R2: numric vector of all observed R2 values and 2) models: a list of the corresponding linear models
calcPhaseShift = function(x,y,left_shift,right_shift)
{
  
  models = list() #list of linear models for each time shift
  R2 = vector(mode = 'numeric') #vector of the correponding R2-values
  
  for(i in left_shift:right_shift)
  {
    if(i<0)
    {
      this.x = x[1:(length(x)-abs(i))]
      this.y = y[(abs(i)+1):length(y)]
    }
    if(i==0)
    {
      this.x = x
      this.y = y
    }
    
    
    if(i>0)
    {
      this.x = x[(abs(i)+1):length(x)]
      this.y = y[1:(length(y)-abs(i))]
    }
    
    mod = lm(this.x~this.y)
    modname = paste('mod',i,sep='_')
    models[[length(models)+1]] <- mod
    #print(summary(mod)) 
    R2 = c(R2,summary(mod)$r.squared)
    
    #rho = cor.test(x,y, method = 'spearman')
    #Rho = c(Rho,rho$estimate)
    
  }
  return(list(R2=R2,models=models))
}
```

## PART I: Field data

This part analyses the field data collected in Königshafen from 2013
until 2019.

Phase shift calculations for Figure 1C,D

``` r
#reducing data set to times of co-occurence
after2015 = subset(M_Mint_prev, M_Mint_prev$Year>2015)

#phase shift
x_phase = -3:4
MintMori_shift = calcPhaseShift(after2015$M_infected_Mint,after2015$M_infected_Mori,-3,4)
bestMod_Index = which(MintMori_shift$R2==max(MintMori_shift$R2))
bestMod = MintMori_shift$models[[bestMod_Index]]
print(paste('Best correlation at shift:',x_phase[bestMod_Index]))
```

    ## [1] "Best correlation at shift: 1"

``` r
summary(bestMod)
```

    ## 
    ## Call:
    ## lm(formula = this.x ~ this.y)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.16786 -0.11268  0.02769  0.10247  0.15117 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.93786    0.08005  11.716 2.57e-06 ***
    ## this.y      -0.71832    0.24064  -2.985   0.0175 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.127 on 8 degrees of freedom
    ## Multiple R-squared:  0.5269, Adjusted R-squared:  0.4678 
    ## F-statistic:  8.91 on 1 and 8 DF,  p-value: 0.01747

Plot Figure 1
B,C,D

![](/README_files/figure-gfm/Plot%201B,C,D-1.png)<!-- -->

Calculate dispersion and reproductive opportunities ‘M\_Mint\_abundance’
is a sampling interval aggregate of M. intestinalis (Mint) abundances in
mussels (M), Mori is M. orientalis, and O is oyster. Different stats are
calculated for each host parasite combination (M\_Mint, M\_mori,
O\_Mori) For all 3 combinations data is merged into M-Mint\_abundance

Analysis for Figure 3B: Aggregation in different parasite-host
combinations

``` r
#analysis
#test for species differences accounting for different variances
fit.gls <- gls(varToMean~hostPar, weights = varIdent(form = ~ 1 | hostPar), data = disp)
summary(fit.gls)
```

    ## Generalized least squares fit by REML
    ##   Model: varToMean ~ hostPar 
    ##   Data: disp 
    ##      AIC      BIC   logLik
    ##   228.77 240.5917 -108.385
    ## 
    ## Variance function:
    ##  Structure: Different standard deviations per stratum
    ##  Formula: ~1 | hostPar 
    ##  Parameter estimates:
    ##   M_Mint   M_Mori   O_Mori 
    ## 1.000000 0.305492 5.467188 
    ## 
    ## Coefficients:
    ##                   Value Std.Error   t-value p-value
    ## (Intercept)    2.971621 0.2644636 11.236411  0.0000
    ## hostParM_Mori -0.980940 0.3056782 -3.209060  0.0023
    ## hostParO_Mori  5.074641 2.7560675  1.841262  0.0712
    ## 
    ##  Correlation: 
    ##               (Intr) hsPM_M
    ## hostParM_Mori -0.865       
    ## hostParO_Mori -0.096  0.083
    ## 
    ## Standardized residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -1.2646412 -0.6527898 -0.3029589  0.3602025  2.7727779 
    ## 
    ## Residual standard error: 1.586781 
    ## Degrees of freedom: 56 total; 53 residual

``` r
anova(fit.gls)
```

    ## Denom. DF: 53 
    ##             numDF   F-value p-value
    ## (Intercept)     1 288.73189  <.0001
    ## hostPar         2   7.38559  0.0015

Analysis Figure 3C: Chances of encountering a mate

``` r
#### prop Repro 2C
#regression before M. orientalis arrival in 2016
M_Mint_abundance_before = subset(M_Mint_abundance,M_Mint_abundance$x_month < 35)
M_Mint_abundance_after = subset(M_Mint_abundance,M_Mint_abundance$x_month > 35)

M_mint_repro_mod_before = lm(M_Mint_abundance_before$repro_M_Mint~ M_Mint_abundance_before$x_month)
print(summary(M_mint_repro_mod_before))
```

    ## 
    ## Call:
    ## lm(formula = M_Mint_abundance_before$repro_M_Mint ~ M_Mint_abundance_before$x_month)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.143510 -0.044232  0.000878  0.057360  0.114619 
    ## 
    ## Coefficients:
    ##                                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                      0.8289533  0.0270430   30.65   <2e-16 ***
    ## M_Mint_abundance_before$x_month -0.0004377  0.0012878   -0.34    0.737    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.07296 on 23 degrees of freedom
    ## Multiple R-squared:  0.004998,   Adjusted R-squared:  -0.03826 
    ## F-statistic: 0.1155 on 1 and 23 DF,  p-value: 0.737

``` r
M_mint_repro_mod_after = lm(M_Mint_abundance_after$repro_M_Mint~ M_Mint_abundance_after$x_month)
print(summary(M_mint_repro_mod_after))
```

    ## 
    ## Call:
    ## lm(formula = M_Mint_abundance_after$repro_M_Mint ~ M_Mint_abundance_after$x_month)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.178387 -0.115694 -0.008255  0.068842  0.308812 
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                     1.118737   0.225927   4.952 0.000789 ***
    ## M_Mint_abundance_after$x_month -0.009723   0.003800  -2.559 0.030751 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1542 on 9 degrees of freedom
    ## Multiple R-squared:  0.4211, Adjusted R-squared:  0.3568 
    ## F-statistic: 6.547 on 1 and 9 DF,  p-value: 0.03075

``` r
M_mori_repro_mod = lm(repro_M_Mori~x_month, data = M_Mint_abundance)
print(summary(M_mori_repro_mod))
```

    ## 
    ## Call:
    ## lm(formula = repro_M_Mori ~ x_month, data = M_Mint_abundance)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.14607 -0.05226  0.02803  0.04729  0.11193 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept) -0.161838   0.148774  -1.088  0.30836   
    ## x_month      0.009051   0.002429   3.726  0.00582 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.08245 on 8 degrees of freedom
    ##   (26 observations deleted due to missingness)
    ## Multiple R-squared:  0.6344, Adjusted R-squared:  0.5887 
    ## F-statistic: 13.88 on 1 and 8 DF,  p-value: 0.00582

``` r
O_mori_repro_mod = lm(repro_O_Mori~x_month, data = M_Mint_abundance)
print(summary(O_mori_repro_mod))
```

    ## 
    ## Call:
    ## lm(formula = repro_O_Mori ~ x_month, data = M_Mint_abundance)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.52545 -0.06460  0.03812  0.12121  0.24152 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept) -0.502173   0.409559  -1.226   0.2550  
    ## x_month      0.018684   0.006687   2.794   0.0234 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.227 on 8 degrees of freedom
    ##   (26 observations deleted due to missingness)
    ## Multiple R-squared:  0.4939, Adjusted R-squared:  0.4306 
    ## F-statistic: 7.807 on 1 and 8 DF,  p-value: 0.02341

Plotting of Figure
3

![](/README_files/figure-gfm/plot%20Figure2-1.png)<!-- -->

## Part II: Simultaneous and Sequential infections (Fig 4)

This part analyses the simultaneous and sequential infection data,
corresponding to Figure 3 \#\#\# Simultaneous infections (Fig
3A)

``` r
#############################################################################################################
# test by binomial GLMM
#############################################################################################################

mod = glmer(cbind(Parasite_count,(exposureDose-Parasite_count))~Parasite * treatment + (1|Mussel_name), data = dataSim, family = 'binomial',control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
print(summary(mod))
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: cbind(Parasite_count, (exposureDose - Parasite_count)) ~ Parasite *  
    ##     treatment + (1 | Mussel_name)
    ##    Data: dataSim
    ## Control: glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+05))
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    256.8    266.3   -123.4    246.8       45 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.1495 -0.2840  0.0019  0.2957  3.1996 
    ## 
    ## Random effects:
    ##  Groups      Name        Variance Std.Dev.
    ##  Mussel_name (Intercept) 0.8514   0.9227  
    ## Number of obs: 50, groups:  Mussel_name, 37
    ## 
    ## Fixed effects:
    ##                           Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                0.04151    0.31002   0.134    0.893    
    ## ParasiteMori              -1.95992    0.30020  -6.529 6.64e-11 ***
    ## treatmentinf              -0.11576    0.42047  -0.275    0.783    
    ## ParasiteMori:treatmentinf  0.07179    0.53450   0.134    0.893    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) PrstMr trtmnt
    ## ParasiteMor -0.333              
    ## treatmentnf -0.737  0.247       
    ## PrstMr:trtm  0.186 -0.545 -0.496

``` r
#############################################################################################################
# pairwise comparisons
#############################################################################################################
dataSim$PxT = interaction(dataSim$Parasite,dataSim$treatment)
compMod = glmer(cbind(Parasite_count,(exposureDose-Parasite_count))~PxT + (1|Mussel_name), data = dataSim, family = 'binomial',control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
print(summary(compMod))
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: cbind(Parasite_count, (exposureDose - Parasite_count)) ~ PxT +  
    ##     (1 | Mussel_name)
    ##    Data: dataSim
    ## Control: glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+05))
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    256.8    266.3   -123.4    246.8       45 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.1495 -0.2840  0.0019  0.2957  3.1996 
    ## 
    ## Random effects:
    ##  Groups      Name        Variance Std.Dev.
    ##  Mussel_name (Intercept) 0.8514   0.9227  
    ## Number of obs: 50, groups:  Mussel_name, 37
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)    0.04151    0.31002   0.134    0.893    
    ## PxTMori.coinf -1.95992    0.30020  -6.529 6.64e-11 ***
    ## PxTMint.inf   -0.11576    0.42047  -0.275    0.783    
    ## PxTMori.inf   -2.00388    0.46620  -4.298 1.72e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) PxTMr.c PxTMn.
    ## PxTMori.cnf -0.333               
    ## PxTMint.inf -0.737  0.247        
    ## PxTMori.inf -0.666  0.242   0.493

``` r
print(summary(glht(compMod, linfct = mcp(PxT = 'Tukey'))))
```

    ## 
    ##   Simultaneous Tests for General Linear Hypotheses
    ## 
    ## Multiple Comparisons of Means: Tukey Contrasts
    ## 
    ## 
    ## Fit: glmer(formula = cbind(Parasite_count, (exposureDose - Parasite_count)) ~ 
    ##     PxT + (1 | Mussel_name), data = dataSim, family = "binomial", 
    ##     control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+05)))
    ## 
    ## Linear Hypotheses:
    ##                              Estimate Std. Error z value Pr(>|z|)    
    ## Mori.coinf - Mint.coinf == 0 -1.95992    0.30020  -6.529  < 1e-04 ***
    ## Mint.inf - Mint.coinf == 0   -0.11576    0.42047  -0.275 0.992431    
    ## Mori.inf - Mint.coinf == 0   -2.00388    0.46620  -4.298 0.000112 ***
    ## Mint.inf - Mori.coinf == 0    1.84416    0.45229   4.077 0.000249 ***
    ## Mori.inf - Mori.coinf == 0   -0.04397    0.48955  -0.090 0.999730    
    ## Mori.inf - Mint.inf == 0     -1.88812    0.44834  -4.211 0.000139 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## (Adjusted p values reported -- single-step method)

``` r
meanInf = aggregate( infRate2 ~ treatment*Parasite, data = dataSim, mean)
meanSE = aggregate( infRate2 ~ treatment*Parasite, data = dataSim, se)
```

### Sequential Infections (Fig 4B)

``` r
#############################################################################################################
# test by binomial GLMM
#############################################################################################################

modSeq = glmer(cbind(inf,(dose-inf))~ Parasite * round * HomHet + (1|Mussel), data = dataSeq, family = 'binomial',control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(modSeq)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: cbind(inf, (dose - inf)) ~ Parasite * round * HomHet + (1 | Mussel)
    ##    Data: dataSeq
    ## Control: glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+05))
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    602.4    638.7   -288.2    576.4      107 
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.53107 -0.58639 -0.05872  0.59034  2.81856 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  Mussel (Intercept) 0.7771   0.8815  
    ## Number of obs: 120, groups:  Mussel, 80
    ## 
    ## Fixed effects:
    ##                                        Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                             0.88115    0.31770   2.774 0.005545 ** 
    ## ParasiteMori                           -1.23849    0.45427  -2.726 0.006405 ** 
    ## roundinf2                              -1.61039    0.46248  -3.482 0.000498 ***
    ## HomHethomolog                          -0.09493    0.47385  -0.200 0.841223    
    ## HomHetheterolog                         0.14151    0.47913   0.295 0.767723    
    ## ParasiteMori:roundinf2                  0.79074    0.65394   1.209 0.226584    
    ## ParasiteMori:HomHethomolog             -1.02098    0.69667  -1.466 0.142780    
    ## ParasiteMori:HomHetheterolog           -1.17788    0.69099  -1.705 0.088263 .  
    ## roundinf2:HomHethomolog                -0.87447    0.57649  -1.517 0.129297    
    ## roundinf2:HomHetheterolog               0.10033    0.68194   0.147 0.883031    
    ## ParasiteMori:roundinf2:HomHethomolog    2.57875    0.82277   3.134 0.001723 ** 
    ## ParasiteMori:roundinf2:HomHetheterolog  1.17617    1.12477   1.046 0.295701    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                      (Intr) PrstMr rndnf2 HmHthm HmHtht PrsM:2 PrstMr:HmHthm
    ## ParasiteMor          -0.699                                                 
    ## roundinf2            -0.689  0.482                                          
    ## HomHethomlg          -0.670  0.468  0.460                                   
    ## HomHethtrlg          -0.662  0.463  0.453  0.444                            
    ## PrstMr:rnd2           0.486 -0.695 -0.705 -0.325 -0.321                     
    ## PrstMr:HmHthm         0.455 -0.651 -0.311 -0.680 -0.303  0.452              
    ## PrstMr:HmHtht         0.458 -0.656 -0.314 -0.308 -0.694  0.456  0.429       
    ## rndnf2:HmHthm         0.551 -0.385 -0.799 -0.542 -0.365  0.565  0.368       
    ## rndnf2:HmHtht         0.466 -0.326 -0.676 -0.312 -0.702  0.477  0.212       
    ## PrstMr:rndnf2:HmHthm -0.385  0.551  0.558  0.380  0.256 -0.794 -0.572       
    ## PrstMr:rndnf2:HmHtht -0.282  0.403  0.409  0.189  0.567 -0.581 -0.263       
    ##                      PrstMr:HmHtht rndnf2:HmHthm rndnf2:HmHtht
    ## ParasiteMor                                                   
    ## roundinf2                                                     
    ## HomHethomlg                                                   
    ## HomHethtrlg                                                   
    ## PrstMr:rnd2                                                   
    ## PrstMr:HmHthm                                                 
    ## PrstMr:HmHtht                                                 
    ## rndnf2:HmHthm         0.253                                   
    ## rndnf2:HmHtht         0.656         0.541                     
    ## PrstMr:rndnf2:HmHthm -0.363        -0.700        -0.379       
    ## PrstMr:rndnf2:HmHtht -0.814        -0.328        -0.809       
    ##                      PrstMr:rndnf2:HmHthm
    ## ParasiteMor                              
    ## roundinf2                                
    ## HomHethomlg                              
    ## HomHethtrlg                              
    ## PrstMr:rnd2                              
    ## PrstMr:HmHthm                            
    ## PrstMr:HmHtht                            
    ## rndnf2:HmHthm                            
    ## rndnf2:HmHtht                            
    ## PrstMr:rndnf2:HmHthm                     
    ## PrstMr:rndnf2:HmHtht  0.462

``` r
#Anova(modSeq, type = 'III')
```

### Pairwise comparisons of Sequential infections

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: cbind(inf, (dose - inf)) ~ PxRxHH + (1 | Mussel)
    ##    Data: dataSeq
    ## Control: glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+05))
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    667.5    698.2   -322.8    645.5      109 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.5314 -0.8444 -0.1775  0.6842  2.8013 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  Mussel (Intercept) 0.7336   0.8565  
    ## Number of obs: 120, groups:  Mussel, 80
    ## 
    ## Fixed effects:
    ##                           Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)               -0.32396    0.30617  -1.058  0.29001   
    ## PxRxHHMint.inf1.heterolog  1.34154    0.46722   2.871  0.00409 **
    ## PxRxHHMint.inf1.single     1.20229    0.43619   2.756  0.00585 **
    ## PxRxHHMint.inf2.heterolog -0.16283    0.46130  -0.353  0.72410   
    ## PxRxHHMint.inf2.single    -0.39906    0.44878  -0.889  0.37389   
    ## PxRxHHMori.homolog        -0.67266    0.44823  -1.501  0.13344   
    ## PxRxHHMori.inf1.heterolog -1.06744    0.48123  -2.218  0.02654 * 
    ## PxRxHHMori.inf1.single    -0.03297    0.44112  -0.075  0.94043   
    ## PxRxHHMori.inf2.heterolog -0.60986    0.47025  -1.297  0.19467   
    ## PxRxHHMori.inf2.single    -0.84732    0.44599  -1.900  0.05745 . 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                 (Intr) PxRxHHMnt.nf1.h PxRxHHMnt.nf1.s PxRxHHMnt.nf2.h
    ## PxRxHHMnt.nf1.h -0.656                                                
    ## PxRxHHMnt.nf1.s -0.702  0.461                                         
    ## PxRxHHMnt.nf2.h -0.664  0.435           0.466                         
    ## PxRxHHMnt.nf2.s -0.682  0.445           0.477           0.453         
    ## PxRxHHMr.hm     -0.683  0.446           0.478           0.453         
    ## PxRxHHMr.nf1.h  -0.636  0.416           0.446           0.761         
    ## PxRxHHMr.nf1.s  -0.694  0.455           0.487           0.461         
    ## PxRxHHMr.nf2.h  -0.651  0.752           0.456           0.432         
    ## PxRxHHMr.nf2.s  -0.686  0.448           0.481           0.456         
    ##                 PxRxHHMnt.nf2.s PxRHHM. PxRxHHMr.nf1.h PxRxHHMr.nf1.s
    ## PxRxHHMnt.nf1.h                                                      
    ## PxRxHHMnt.nf1.s                                                      
    ## PxRxHHMnt.nf2.h                                                      
    ## PxRxHHMnt.nf2.s                                                      
    ## PxRxHHMr.hm      0.468                                               
    ## PxRxHHMr.nf1.h   0.435           0.435                               
    ## PxRxHHMr.nf1.s   0.474           0.474   0.442                       
    ## PxRxHHMr.nf2.h   0.445           0.446   0.415          0.452        
    ## PxRxHHMr.nf2.s   0.471           0.471   0.437          0.476        
    ##                 PxRxHHMr.nf2.h
    ## PxRxHHMnt.nf1.h               
    ## PxRxHHMnt.nf1.s               
    ## PxRxHHMnt.nf2.h               
    ## PxRxHHMnt.nf2.s               
    ## PxRxHHMr.hm                   
    ## PxRxHHMr.nf1.h                
    ## PxRxHHMr.nf1.s                
    ## PxRxHHMr.nf2.h                
    ## PxRxHHMr.nf2.s   0.448

    ## Warning in RET$pfunction("adjusted", ...): Completion with error > abseps
    
    ## Warning in RET$pfunction("adjusted", ...): Completion with error > abseps
    
    ## Warning in RET$pfunction("adjusted", ...): Completion with error > abseps
    
    ## Warning in RET$pfunction("adjusted", ...): Completion with error > abseps
    
    ## Warning in RET$pfunction("adjusted", ...): Completion with error > abseps

    ## 
    ##   Simultaneous Tests for General Linear Hypotheses
    ## 
    ## Multiple Comparisons of Means: Tukey Contrasts
    ## 
    ## 
    ## Fit: glmer(formula = cbind(inf, (dose - inf)) ~ PxRxHH + (1 | Mussel), 
    ##     data = dataSeq, family = "binomial", control = glmerControl(optimizer = "bobyqa", 
    ##         optCtrl = list(maxfun = 2e+05)))
    ## 
    ## Linear Hypotheses:
    ##                                                Estimate Std. Error z value
    ## Mint.inf1.heterolog - Mint.homolog == 0         1.34154    0.46722   2.871
    ## Mint.inf1.single - Mint.homolog == 0            1.20229    0.43619   2.756
    ## Mint.inf2.heterolog - Mint.homolog == 0        -0.16283    0.46130  -0.353
    ## Mint.inf2.single - Mint.homolog == 0           -0.39906    0.44878  -0.889
    ## Mori.homolog - Mint.homolog == 0               -0.67266    0.44823  -1.501
    ## Mori.inf1.heterolog - Mint.homolog == 0        -1.06744    0.48123  -2.218
    ## Mori.inf1.single - Mint.homolog == 0           -0.03297    0.44112  -0.075
    ## Mori.inf2.heterolog - Mint.homolog == 0        -0.60986    0.47025  -1.297
    ## Mori.inf2.single - Mint.homolog == 0           -0.84732    0.44599  -1.900
    ## Mint.inf1.single - Mint.inf1.heterolog == 0    -0.13925    0.46959  -0.297
    ## Mint.inf2.heterolog - Mint.inf1.heterolog == 0 -1.50437    0.49359  -3.048
    ## Mint.inf2.single - Mint.inf1.heterolog == 0    -1.74060    0.48293  -3.604
    ## Mori.homolog - Mint.inf1.heterolog == 0        -2.01420    0.48217  -4.177
    ## Mori.inf1.heterolog - Mint.inf1.heterolog == 0 -2.40898    0.51256  -4.700
    ## Mori.inf1.single - Mint.inf1.heterolog == 0    -1.37451    0.47476  -2.895
    ## Mori.inf2.heterolog - Mint.inf1.heterolog == 0 -1.95140    0.33020  -5.910
    ## Mori.inf2.single - Mint.inf1.heterolog == 0    -2.18886    0.48024  -4.558
    ## Mint.inf2.heterolog - Mint.inf1.single == 0    -1.36512    0.46431  -2.940
    ## Mint.inf2.single - Mint.inf1.single == 0       -1.60135    0.45250  -3.539
    ## Mori.homolog - Mint.inf1.single == 0           -1.87495    0.45180  -4.150
    ## Mori.inf1.heterolog - Mint.inf1.single == 0    -2.26973    0.48430  -4.687
    ## Mori.inf1.single - Mint.inf1.single == 0       -1.23526    0.44426  -2.780
    ## Mori.inf2.heterolog - Mint.inf1.single == 0    -1.81215    0.47348  -3.827
    ## Mori.inf2.single - Mint.inf1.single == 0       -2.04961    0.44967  -4.558
    ## Mint.inf2.single - Mint.inf2.heterolog == 0    -0.23623    0.47613  -0.496
    ## Mori.homolog - Mint.inf2.heterolog == 0        -0.50983    0.47562  -1.072
    ## Mori.inf1.heterolog - Mint.inf2.heterolog == 0 -0.90461    0.32616  -2.774
    ## Mori.inf1.single - Mint.inf2.heterolog == 0     0.12986    0.46893   0.277
    ## Mori.inf2.heterolog - Mint.inf2.heterolog == 0 -0.44703    0.49643  -0.900
    ## Mori.inf2.single - Mint.inf2.heterolog == 0    -0.68449    0.47350  -1.446
    ## Mori.homolog - Mint.inf2.single == 0           -0.27359    0.46268  -0.591
    ## Mori.inf1.heterolog - Mint.inf2.single == 0    -0.66838    0.49516  -1.350
    ## Mori.inf1.single - Mint.inf2.single == 0        0.36609    0.45663   0.802
    ## Mori.inf2.heterolog - Mint.inf2.single == 0    -0.21080    0.48434  -0.435
    ## Mori.inf2.single - Mint.inf2.single == 0       -0.44826    0.46034  -0.974
    ## Mori.inf1.heterolog - Mori.homolog == 0        -0.39478    0.49474  -0.798
    ## Mori.inf1.single - Mori.homolog == 0            0.63969    0.45609   1.403
    ## Mori.inf2.heterolog - Mori.homolog == 0         0.06280    0.48394   0.130
    ## Mori.inf2.single - Mori.homolog == 0           -0.17466    0.46004  -0.380
    ## Mori.inf1.single - Mori.inf1.heterolog == 0     1.03447    0.48856   2.117
    ## Mori.inf2.heterolog - Mori.inf1.heterolog == 0  0.45758    0.51486   0.889
    ## Mori.inf2.single - Mori.inf1.heterolog == 0     0.22012    0.49266   0.447
    ## Mori.inf2.heterolog - Mori.inf1.single == 0    -0.57689    0.47774  -1.208
    ## Mori.inf2.single - Mori.inf1.single == 0       -0.81435    0.45388  -1.794
    ## Mori.inf2.single - Mori.inf2.heterolog == 0    -0.23746    0.48179  -0.493
    ##                                                Pr(>|z|)    
    ## Mint.inf1.heterolog - Mint.homolog == 0          0.1091    
    ## Mint.inf1.single - Mint.homolog == 0             0.1461    
    ## Mint.inf2.heterolog - Mint.homolog == 0          1.0000    
    ## Mint.inf2.single - Mint.homolog == 0             0.9965    
    ## Mori.homolog - Mint.homolog == 0                 0.8860    
    ## Mori.inf1.heterolog - Mint.homolog == 0          0.4317    
    ## Mori.inf1.single - Mint.homolog == 0             1.0000    
    ## Mori.inf2.heterolog - Mint.homolog == 0          0.9514    
    ## Mori.inf2.single - Mint.homolog == 0             0.6573    
    ## Mint.inf1.single - Mint.inf1.heterolog == 0      1.0000    
    ## Mint.inf2.heterolog - Mint.inf1.heterolog == 0   0.0677 .  
    ## Mint.inf2.single - Mint.inf1.heterolog == 0      0.0111 *  
    ## Mori.homolog - Mint.inf1.heterolog == 0           <0.01 ** 
    ## Mori.inf1.heterolog - Mint.inf1.heterolog == 0    <0.01 ***
    ## Mori.inf1.single - Mint.inf1.heterolog == 0      0.1020    
    ## Mori.inf2.heterolog - Mint.inf1.heterolog == 0    <0.01 ***
    ## Mori.inf2.single - Mint.inf1.heterolog == 0       <0.01 ***
    ## Mint.inf2.heterolog - Mint.inf1.single == 0      0.0912 .  
    ## Mint.inf2.single - Mint.inf1.single == 0         0.0142 *  
    ## Mori.homolog - Mint.inf1.single == 0              <0.01 ** 
    ## Mori.inf1.heterolog - Mint.inf1.single == 0       <0.01 ***
    ## Mori.inf1.single - Mint.inf1.single == 0         0.1375    
    ## Mori.inf2.heterolog - Mint.inf1.single == 0       <0.01 ** 
    ## Mori.inf2.single - Mint.inf1.single == 0          <0.01 ***
    ## Mint.inf2.single - Mint.inf2.heterolog == 0      1.0000    
    ## Mori.homolog - Mint.inf2.heterolog == 0          0.9863    
    ## Mori.inf1.heterolog - Mint.inf2.heterolog == 0   0.1403    
    ## Mori.inf1.single - Mint.inf2.heterolog == 0      1.0000    
    ## Mori.inf2.heterolog - Mint.inf2.heterolog == 0   0.9962    
    ## Mori.inf2.single - Mint.inf2.heterolog == 0      0.9074    
    ## Mori.homolog - Mint.inf2.single == 0             0.9999    
    ## Mori.inf1.heterolog - Mint.inf2.single == 0      0.9378    
    ## Mori.inf1.single - Mint.inf2.single == 0         0.9984    
    ## Mori.inf2.heterolog - Mint.inf2.single == 0      1.0000    
    ## Mori.inf2.single - Mint.inf2.single == 0         0.9931    
    ## Mori.inf1.heterolog - Mori.homolog == 0          0.9985    
    ## Mori.inf1.single - Mori.homolog == 0             0.9221    
    ## Mori.inf2.heterolog - Mori.homolog == 0          1.0000    
    ## Mori.inf2.single - Mori.homolog == 0             1.0000    
    ## Mori.inf1.single - Mori.inf1.heterolog == 0      0.5025    
    ## Mori.inf2.heterolog - Mori.inf1.heterolog == 0   0.9965    
    ## Mori.inf2.single - Mori.inf1.heterolog == 0      1.0000    
    ## Mori.inf2.heterolog - Mori.inf1.single == 0      0.9692    
    ## Mori.inf2.single - Mori.inf1.single == 0         0.7286    
    ## Mori.inf2.single - Mori.inf2.heterolog == 0      1.0000    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## (Adjusted p values reported -- single-step method)

Plotting Figure
    4

    ## [1] "Using c(0.5,4) as input for xlim, note that default values for these dimensions are c(0.5,2.5)"

    ## [1] "Using c(0.5,10.5) as input for xlim, note that default values for these dimensions are c(0.5,7)"

![](/README_files/figure-gfm/Figure3-1.png)<!-- -->

## PART III: Host dependence of infectivity

Tests whether the source species (mussel or oyster) has an effect on the
infectitivity of parasite offspring in both species

Testing for cost of
generalism

``` r
parasite.mod = glm(cbind(M_Total,(Cop_infiltrated-M_Total))~HostParasite, data = dataHostNoMiCG, family = 'binomial')
summary(parasite.mod)
```

    ## 
    ## Call:
    ## glm(formula = cbind(M_Total, (Cop_infiltrated - M_Total)) ~ HostParasite, 
    ##     family = "binomial", data = dataHostNoMiCG)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.2296  -2.2090  -0.3385   1.0575   5.0613  
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)               0.4017     0.1369   2.933  0.00335 ** 
    ## HostParasiteMo.M_edulis  -2.6795     0.2031 -13.192  < 2e-16 ***
    ## HostParasiteMo.C_gigas   -2.5579     0.2061 -12.409  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 456.63  on 50  degrees of freedom
    ## Residual deviance: 212.94  on 48  degrees of freedom
    ## AIC: 316.49
    ## 
    ## Number of Fisher Scoring iterations: 5

Testing for host source
dependence

``` r
parental.mod = glm(cbind(M_Total,(Cop_infiltrated-M_Total))~Treatment, data = dataHostNoMiCG, family = 'binomial')
summary(parental.mod)
```

    ## 
    ## Call:
    ## glm(formula = cbind(M_Total, (Cop_infiltrated - M_Total)) ~ Treatment, 
    ##     family = "binomial", data = dataHostNoMiCG)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.2296  -1.6224  -0.8825   0.8148   5.0613  
    ## 
    ## Coefficients:
    ##               Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)     0.4017     0.1369   2.933  0.00335 ** 
    ## TreatmentMoCg  -2.1522     0.1885 -11.419  < 2e-16 ***
    ## TreatmentMoMe  -3.3195     0.2435 -13.632  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 456.63  on 50  degrees of freedom
    ## Residual deviance: 186.61  on 48  degrees of freedom
    ## AIC: 290.16
    ## 
    ## Number of Fisher Scoring iterations: 5

``` r
anova(parental.mod,test = 'LRT')
```

    ## Analysis of Deviance Table
    ## 
    ## Model: binomial, link: logit
    ## 
    ## Response: cbind(M_Total, (Cop_infiltrated - M_Total))
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ## 
    ##           Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
    ## NULL                         50     456.63              
    ## Treatment  2   270.02        48     186.61 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
dataHost_Mori = subset(dataHostNoMiCG,dataHostNoMiCG$Treatment != 'Mi')
hostMori.mod = glm(cbind(M_Total,(Cop_infiltrated-M_Total))~Treatment * HostSpecies, data = dataHost_Mori, family = 'binomial')
summary(hostMori.mod)
```

    ## 
    ## Call:
    ## glm(formula = cbind(M_Total, (Cop_infiltrated - M_Total)) ~ Treatment * 
    ##     HostSpecies, family = "binomial", data = dataHost_Mori)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -2.8652  -1.6343  -0.9588   0.8979   3.7891  
    ## 
    ## Coefficients:
    ##                                   Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                        -1.6792     0.1815  -9.249  < 2e-16 ***
    ## TreatmentMoMe                      -1.2842     0.3586  -3.582 0.000341 ***
    ## HostSpeciesM_edulis                -0.1416     0.2591  -0.546 0.584782    
    ## TreatmentMoMe:HostSpeciesM_edulis   0.2220     0.4828   0.460 0.645620    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 138.63  on 41  degrees of freedom
    ## Residual deviance: 111.64  on 38  degrees of freedom
    ## AIC: 192.06
    ## 
    ## Number of Fisher Scoring iterations: 5

Plotting Figure
5

![](/README_files/figure-gfm/Figure%204-1.png)<!-- -->

## PART IV: Infection preference

tests whether the prior presence of other con- or heterospecific
infections changes the infection choice when offered an infected vs. an
uninfected
host

``` r
# linear model with infection choice parasite species and homologous vs heterologous infection as factor
mod_hom_het = lm(SecondInfection_prop_difference~ FirstIntensity + SecondInfection * hom_het  , data = mytilus_noPreinf)
#anova table
anova(mod_hom_het)
```

    ## Analysis of Variance Table
    ## 
    ## Response: SecondInfection_prop_difference
    ##                         Df Sum Sq Mean Sq F value   Pr(>F)   
    ## FirstIntensity           1 0.0061 0.00608  0.0720 0.789893   
    ## SecondInfection          1 0.8918 0.89177 10.5655 0.002416 **
    ## hom_het                  1 0.0198 0.01979  0.2345 0.631022   
    ## SecondInfection:hom_het  1 0.4427 0.44268  5.2448 0.027654 * 
    ## Residuals               38 3.2073 0.08440                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(mod_hom_het)
```

    ## 
    ## Call:
    ## lm(formula = SecondInfection_prop_difference ~ FirstIntensity + 
    ##     SecondInfection * hom_het, data = mytilus_noPreinf)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.72467 -0.17073 -0.03532  0.15615  0.65737 
    ## 
    ## Coefficients:
    ##                                             Estimate Std. Error t value
    ## (Intercept)                                 -0.15586    0.09608  -1.622
    ## FirstIntensity                              -0.01690    0.01082  -1.562
    ## SecondInfectionOrientalis                    0.53954    0.14015   3.850
    ## hom_hethomologous                            0.28950    0.13981   2.071
    ## SecondInfectionOrientalis:hom_hethomologous -0.47510    0.20745  -2.290
    ##                                             Pr(>|t|)    
    ## (Intercept)                                  0.11304    
    ## FirstIntensity                               0.12659    
    ## SecondInfectionOrientalis                    0.00044 ***
    ## hom_hethomologous                            0.04524 *  
    ## SecondInfectionOrientalis:hom_hethomologous  0.02765 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2905 on 38 degrees of freedom
    ## Multiple R-squared:  0.2978, Adjusted R-squared:  0.2239 
    ## F-statistic: 4.029 on 4 and 38 DF,  p-value: 0.008053

Plot Figure
6

![](/README_files/figure-gfm/Figure%205-1.png)<!-- -->
