library(raster)

## yogyakarta location
yogyakarta.x <- 110.36444
yogyakarta.y <- -7.801389

## get raster of global FOI from Cattarino2020 paper
foi.raster <- raster("~/Documents/awed_trial_modeling/aax4144_data_file_s3.tif")
foi.df <- as.data.frame(foi.raster,xy=T)
coord.x <- foi.df$x[which.min(abs(foi.df$x - yogyakarta.x))]
coord.y <- foi.df$y[which.min(abs(foi.df$y - yogyakarta.y))]
foi.val <- unlist(foi.df[foi.df$x == coord.x & foi.df$y == coord.y,3])

## get first R0 raster
R0a.raster <- raster("~/Documents/awed_trial_modeling/aax4144_data_file_s4.tif")
R0a.df <- as.data.frame(R0a.raster,xy=T)
coord.x <- R0a.df$x[which.min(abs(R0a.df$x - yogyakarta.x))]
coord.y <- R0a.df$y[which.min(abs(R0a.df$y - yogyakarta.y))]
R0a.val <- unlist(R0a.df[R0a.df$x == coord.x & R0a.df$y == coord.y,3])

## get second R0 raster
R0b.raster <- raster("~/Documents/awed_trial_modeling/aax4144_data_file_s5.tif")
R0b.df <- as.data.frame(R0b.raster,xy=T)
coord.x <- R0b.df$x[which.min(abs(R0b.df$x - yogyakarta.x))]
coord.y <- R0b.df$y[which.min(abs(R0b.df$y - yogyakarta.y))]
R0b.val <- unlist(R0b.df[R0b.df$x == coord.x & R0b.df$y == coord.y,3])

## save
save(foi.val,R0a.val,R0b.val,file="./foi_and_r0.RData")
