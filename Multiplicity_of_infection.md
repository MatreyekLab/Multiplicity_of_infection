Multiplicity\_of\_infection
================
Kenneth Matreyek
8/22/2020

``` r
library(tidyverse)
```

    ## ── Attaching packages ────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.2     ✓ purrr   0.3.4
    ## ✓ tibble  3.0.1     ✓ dplyr   1.0.0
    ## ✓ tidyr   1.1.0     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.5.0

    ## ── Conflicts ───────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
test_fxn <- function(lambda,k)(lambda^k*exp(-lambda))/factorial(k) #lambda is moi, and k is singly, doubly, etc infected.
test_fxn(1,2)
```

    ## [1] 0.1839397

``` r
lambda_vector <- c(seq(0.01,0.1,0.01),seq(0.1,1,0.1),seq(1,10,1))
instance_type <- seq(0,20)

test_frame <- data.frame(matrix(nrow = length(lambda_vector), ncol = length(instance_type),0))
#rownames(test_frame) <- paste("moi_",lambda_vector,sep="")
colnames(test_frame) <- instance_type

for(lambda in 1:nrow(test_frame)){
  for(k in 1:ncol(test_frame)){
    test_frame[lambda,k] <- round(test_fxn(lambda_vector[lambda],instance_type[k]),3)
  }
}

summary_frame <- data.frame(cbind(test_frame[,1],rowSums(test_frame[,2:ncol(test_frame)])))
colnames(summary_frame) <- c("uninfected","infected")
summary_frame$lambda <- lambda_vector

MOI_Plot <- ggplot() + theme_bw() + scale_x_log10() + scale_y_log10(limits = c(0.1,100)) +
  labs(x = "Multiplicity of infection", y = "Percent of infected cells") +
  geom_hline(yintercept = 100, size = 2, alpha = 0.4) +
  geom_hline(yintercept = 20, size = 1, linetype = 2) +
  geom_line(data = summary_frame, aes(x = lambda, y = infected*100), alpha = 0.4) +
  geom_point(data = summary_frame, aes(x = lambda, y = infected*100))
MOI_Plot
```

    ## Warning: Removed 1 rows containing missing values (geom_point).

![](Multiplicity_of_infection_files/figure-gfm/MOI%20analysis-1.png)<!-- -->

``` r
ggsave(file = "MOI_Plot2.png", MOI_Plot, height = 4, width = 4)
```

    ## Warning: Removed 1 rows containing missing values (geom_point).

``` r
# Pretending we're running 100,000 cells in the flow

sampling_frame <- data.frame("prob" = c(seq(0.0005,0.001,0.0001),seq(0.001,0.01,0.001),
                                        seq(0.01,0.1,0.01),seq(0.1,1,0.1),seq(1,10,1),20))
sampling_frame$cv <- 0
sampling_frame$pct_uninform <- 0

for(x in 1:nrow(sampling_frame)){
  temp_frame <- rbinom(size = 300000, prob = sampling_frame$prob[x]/100, n = 100)
  sampling_frame$cv[x] <- sd(temp_frame, na.rm = T)/mean(temp_frame, na.rm = T)
  if(is.na(sampling_frame$cv[x])){sampling_frame$cv[x] <- 0}
  sampling_frame$pct_uninform[x] <- sum(temp_frame == 0)/length(temp_frame)
}

sampling_frame$pass_criteria <- "no"
for(x in 1:nrow(sampling_frame)){if(sampling_frame$cv[x] <= 0.5 & sampling_frame$pct_uninform[x] <= 0.05){sampling_frame$pass_criteria[x] <- "yes";break}}


CV_plot <- ggplot() + theme_bw() + scale_x_log10() + scale_y_log10() +
  labs(x = "Percent GFP positive", y = "Coefficient of Variaiton") +
  geom_line(data = sampling_frame, aes(x = prob, y = cv), alpha = 0.4) +
  geom_point(data = sampling_frame, aes(x = prob, y = cv)) +
  geom_vline(xintercept = subset(sampling_frame, pass_criteria == "yes")[,"prob"], linetype = 2)
CV_plot
```

![](Multiplicity_of_infection_files/figure-gfm/MOI%20analysis-2.png)<!-- -->

``` r
#ggsave(file = "MOI_Plot2.png", MOI_Plot, height = 4, width = 4)

uninform_plot <- ggplot() + theme_bw() + scale_x_log10() + scale_y_log10() +
  labs(x = "Percent GFP positive", y = "Fraction of replicates with zero counts") +
  geom_line(data = sampling_frame, aes(x = prob, y = pct_uninform), alpha = 0.4) +
  geom_point(data = sampling_frame, aes(x = prob, y = pct_uninform)) +
  geom_vline(xintercept = subset(sampling_frame, pass_criteria == "yes")[,"prob"], linetype = 2)
uninform_plot
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](Multiplicity_of_infection_files/figure-gfm/MOI%20analysis-3.png)<!-- -->

``` r
#ggsave(file = "Uninform_plot", uninform_plot, height = 4, width = 4)

library(patchwork)
combined_plot <- CV_plot / uninform_plot
combined_plot
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    
    ## Warning: Transformation introduced infinite values in continuous y-axis

![](Multiplicity_of_infection_files/figure-gfm/MOI%20analysis-4.png)<!-- -->

``` r
ggsave(file = "Combined_plot.png", combined_plot, height = 4*1.25, width = 6*1.25)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    
    ## Warning: Transformation introduced infinite values in continuous y-axis
