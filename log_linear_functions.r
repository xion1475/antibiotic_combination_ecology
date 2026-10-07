library(zoo)
library(ggplot2)
library(dplyr)
library("ggpubr")
library(scales)
library("openxlsx")
library(ggsignif)
library(tidyverse)
library(patchwork)

fit_loglinear = function(x, y, tries = 5){
  m1 <- NA
  #figure out start lag guess
  ymax = max(y)
  ymin = min(y)
  saturation_time = x[y > (ymin + ((ymax - ymin) * 0.95))][1]
  lag_time = x[y > (ymin + ((ymax - ymin) * 0.1))][1]
  while(is.na(m1[1]) && (tries > 0)){
    m1 <- tryCatch( nls(y ~ log_linear(x, r, y0, lag, end),
                        start = list(r = runif(1, 0.1, 2), 
                                     y0 = min(y), 
                                     lag = lag_time, 
                                     end = saturation_time)),
                    error = function(e) return(NA))
    tries = tries - 1
  }
  if (is.na(m1[1])){
    return(c(NA, NA,NA,NA))
  }
  r = coef(m1)[1]
  y0 = coef(m1)[2]
  start = coef(m1)[3]
  end = coef(m1)[4]
  return(c(r, y0, start, end))
}


fit_logistic = function(x, y, tries = 5){
    m1 <- NA
    #figure out start lag guess
    lag_guess = guess_half_max(x,y)
    while(is.na(m1[1]) && (tries > 0)){
      m1 <- tryCatch( nls(y ~ logistic(x, r, lag, K, y0),
                          start = list(r = runif(1, 0.1, 1), 
                                       lag = lag_guess + sample(-20:20, 1), 
                                       K = max(y), 
                                       y0 = min(y))),
                      error = function(e) return(NA))
      tries = tries - 1
    }
    if (is.na(m1)){
      return(c(NA, NA,NA,NA))
    }
    #   m1 <- nls(y ~ baranyi(x, r, lag, ymax, y0),
    #             start = list(r = 0.1, lag = x[round(length(x) / 3)], ymax = max(y), y0 = min(y)))
    r = coef(m1)[1]
    lag = coef(m1)[2]
    K = coef(m1)[3]
    y0 = coef(m1)[4]
    return(c(r, lag, K, y0))
}

fit_baranyi = function(x, y, tries = 100){
  m1 <- NA
  #figure out start lag guess
  lag_guess = guess_half_max(x,y)
  while(is.na(m1[1]) && (tries > 0)){
    m1 <- tryCatch( nls(y ~ baranyi(x, r, lag, ymax, y0),
                        start = list(r = runif(1, 0.1, 0.4), 
                                     lag = lag_guess + sample(-20:20, 1), 
                                     ymax = max(y), 
                                     y0 = min(y))),
                    error = function(e) return(NA))
    tries = tries - 1
  }
  if (is.na(m1[1])){
    return(c(NA, NA,NA,NA))
  }
  #   m1 <- nls(y ~ baranyi(x, r, lag, ymax, y0),
  #             start = list(r = 0.1, lag = x[round(length(x) / 3)], ymax = max(y), y0 = min(y)))
  r = coef(m1)[1]
  lag = coef(m1)[2]
  ymax = coef(m1)[3]
  y0 = coef(m1)[4]
  return(c(r, lag, ymax, y0))
}

guess_half_max = function(x, y){
  # this function looks at a logistically-increasing time series
  # and guesses when the growth rate is at a maximum (which is also when
  # the trend is at half-max)
  
  # put in order
  y = y[sort(x, index.return = TRUE)$ix]
  x = sort(x)
  
  #find approximate time of max growth rate, using diffs, on smoothed y
  y2 = diff(y)
  y2 = zoo::rollmean(y2, 11, fill = NA, align = "center")
  half_max_idx = which(y2 == max(y2, na.rm = TRUE))[1]
  half_max = x[half_max_idx]
  return(half_max)
}

baranyi <- function(t, r, lag, logymax, logy0){
  At = t + (1 / r) * log(exp(-r * t) + exp(-r * lag) - exp(-r * (t + lag)))
  logy = logy0 + r * At - log(1 + (exp(r * At) - 1) / exp(logymax - logy0))
  return(logy)
}

logistic <- function(t, r, lag, K, y0){
  y0 + (K - y0) / (1 + exp(-r * (t-lag)))
}

log_linear <- function(t, r, y0, lag, end){
  y = rep(0, length(t))
  y[t < lag] = y0
  y[t >= lag & t < end] = y0 + (t[t >= lag & t < end]-lag) * r
  y[t >= end] = y0 + (end - lag) * r
  return(y)
}

data_normalization_by_first_data_point = function (OD_path, cutoff) {
  raw_data=read.delim(OD_path, header = FALSE, sep = "\t", row.names = 1)

  finalized_data=raw_data
  total_time_point=ncol(raw_data)
  for (i in 4:nrow(raw_data)) {
    row_name = rownames(finalized_data[i,])
    finalized_data[i,] = raw_data[i,] - raw_data[i,1]
    finalized_data[i,] = ifelse(finalized_data[i,] < 0, cutoff, finalized_data[i,] + cutoff)
    #finalized_data[i,] = ifelse(finalized_data[i,] < cutoff, cutoff, finalized_data[i,])
  }

  return(finalized_data)
}

data_normalization_by_max_average_rate = function(finalized_data, turn_type) {
  
  total_time_point=ncol(finalized_data)
  for (i in 4:nrow(finalized_data)) {
    rate_list = c(NA, diff(t(finalized_data[i,])))  
    max_idx <- which.max(rate_list)[1]
    
    if (turn_type == "negative") {
      turn_idx <- which(rate_list <= 0 & seq_along(rate_list) > max_idx)[1]
    } else if (turn_type == "positive") {
      threshold <- mean(rate_list, na.rm = TRUE)
      turn_idx <- which(rate_list < threshold & seq_along(rate_list) > max_idx)[1]
    }
    finalized_data[i,turn_idx:total_time_point] = finalized_data[i, turn_idx-1]
  }
  return(finalized_data)
}
