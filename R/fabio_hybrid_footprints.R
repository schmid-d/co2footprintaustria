##############################################################################################
##  FABIO hybrid footprints
##############################################################################################

library(Matrix)
library(tidyverse)

rm(list=ls()); gc()

is.finite.data.frame <- function(x) do.call(cbind, lapply(x, is.finite))
agg <- function(x) { x <- as.matrix(x) %*% sapply(unique(colnames(x)),"==",colnames(x));  return(x) }


#-------------------------------------------------------------------------
# Make intitial settings
#-------------------------------------------------------------------------
# read region classification
regions <- read.csv(file="./input/fabio_countries.csv", header=TRUE, stringsAsFactors = FALSE)
# read commodity classification
items <- read.csv(file="./input/fabio_items.csv", header=TRUE)
nrreg <- nrow(regions)
nrcom <- nrow(items)
index <- data.frame(ISO = rep(regions$ISO, each = nrcom),
                    country = rep(regions$Country, each = nrcom),
                    item = rep(items$Item, nrreg),
                    group = rep(items$Group, nrreg))

allocation <- c("mass","value")[1]
year <- 2012
#-------------------------------------------------------------------------
# Read data
#-------------------------------------------------------------------------
if(allocation=="mass") L <- readRDS(file=paste0("/mnt/nfs_fineprint/tmp/fabio/hybrid/",year,"_B_inv_mass.rds"))
if(allocation=="value") L <- readRDS(file=paste0("/mnt/nfs_fineprint/tmp/fabio/hybrid/",year,"_B_inv_price.rds"))

X <- readRDS(file=paste0("/mnt/nfs_fineprint/tmp/fabio/",year,"_X.rds"))
E <- readRDS(file=paste0("/mnt/nfs_fineprint/tmp/fabio/",year,"_E.rds"))

grazing <- read.csv("input/grazing.csv")
E$Landuse[E$Item.Code==2001] <- E$Biomass[E$Item.Code==2001] / grazing$t_per_ha

load(file=paste0("/mnt/nfs_fineprint/tmp/exiobase/pxp/",year,"_x.RData"))
load(file=paste0("/mnt/nfs_fineprint/tmp/exiobase/pxp/",year,"_Y.RData"))

load(file="/mnt/nfs_fineprint/tmp/exiobase/Y.codes.RData")
load(file="/mnt/nfs_fineprint/tmp/exiobase/pxp/IO.codes.RData")


footprint <- function(country = "EU", extension = "Landuse", allocation = "value"){
  #-------------------------------------------------------------------------
  # Prepare Multipliers
  #-------------------------------------------------------------------------
  ext <- as.vector(E[,extension]) / X
  ext[!is.finite(ext)] <- 0
  ext[ext < 0] <- 0         # eliminate negative values
  MP <- ext * L
  #-------------------------------------------------------------------------
  # Prepare Final Demand
  #-------------------------------------------------------------------------
  if(country=="EU27"){
    Y_country <- Y[,Y.codes$`Region Name` %in% unique(Y.codes$`Region Name`)[1:27]]
    colnames(Y_country) <- Y.codes$`Final Demand Category`[Y.codes$`Region Name` %in% unique(Y.codes$`Region Name`)[1:27]]
    Y_country <- agg(Y_country)
  } else if(country=="EU"){
    Y_country <- Y[,Y.codes$`Region Name` %in% unique(Y.codes$`Region Name`)[1:28]]
    colnames(Y_country) <- Y.codes$`Final Demand Category`[Y.codes$`Region Name` %in% unique(Y.codes$`Region Name`)[1:28]]
    Y_country <- agg(Y_country)
  } else {
    Y_country <- Y[,Y.codes$`Region Name` == country]
    colnames(Y_country) <- Y.codes$`Final Demand Category`[Y.codes$`Region Name` == country]
  }
  #-------------------------------------------------------------------------
  # Calculate detailed Footprints
  #-------------------------------------------------------------------------
  FP <- t(t(MP) * rowSums(Y_country))
  colnames(FP) <- paste0(IO.codes$Country.Code, "_", IO.codes$Product.Name)
  rownames(FP) <- paste0(index$ISO, "_", index$item)
  results <- FP %>%
    as.matrix() %>% 
    as_tibble() %>% 
    mutate(origin = paste0(index$ISO, "_", index$item)) %>% 
    gather(index, value, -origin) %>%
    mutate(country_origin = substr(origin,1,3)) %>% 
    mutate(item_origin = substr(origin,5,100)) %>% 
    mutate(country_target = substr(index,1,2)) %>% 
    mutate(item_target = substr(index,4,100)) %>% 
    select(-index, -origin) %>% 
    filter(value != 0)
  
  results$group_origin <- items$Com.Group[match(results$item_origin,items$Item)]
  results$group_target <- IO.codes$Sector.Group[match(results$item_target, IO.codes$Product.Name)]
  
  # data.table::fwrite(results, file=paste0("./output/FABIO_hybrid_",country,"_",year,"_",extension,"_results_",allocation,"-alloc_full.csv"), sep=",")
  
  data <- results %>% 
    group_by(item_target, group_origin) %>% 
    summarise(value = round(sum(value))) %>% 
    filter(value != 0) %>% 
    spread(group_origin, value, fill = 0)
  data$index <- IO.codes$Index[IO.codes$Country.Code=="AT"][match(substr(data$item_target,1,50), substr(IO.codes$Product.Name[IO.codes$Country.Code=="AT"],1,50))]
  data <- data[,c(ncol(data),1:(ncol(data)-1))]
  data <- data[order(data$index),]
  data.table::fwrite(data, file=paste0("./output/FABIO-hybrid_",country,"_",year,"_",extension,"_OtherUses_",allocation,"-alloc.csv"), sep=",")
  
  results$region_origin <- regions$EU27[match(results$country_origin, regions$ISO)]
  if(! country %in% regions$EU27){
    results$region_origin[results$country_origin == regions$ISO[regions$ISO2==country & !is.na(regions$ISO2)]] <- 
      regions$ISO[regions$ISO2==country & !is.na(regions$ISO2)]
  }
  
  data <- results %>% 
    group_by(item_target, region_origin) %>% 
    summarise(value = round(sum(value))) %>% 
    filter(value != 0) %>% 
    spread(region_origin, value, fill = 0)
  data.table::fwrite(data, file=paste0("./output/FABIO-hybrid_",country,"_",year,"_",extension,"_OtherUses_",allocation,"-alloc_origin-final-product.csv"), sep=",")
  
  data <- results %>% 
    group_by(item_origin, region_origin) %>% 
    summarise(value = round(sum(value))) %>% 
    filter(value != 0) %>% 
    spread(region_origin, value, fill = 0)
  data.table::fwrite(data, file=paste0("./output/FABIO-hybrid_",country,"_",year,"_",extension,"_OtherUses_",allocation,"-alloc_origin-primary-crops.csv"), sep=",")
  
  return(data)
}


#-------------------------------------------------------------------------
# Calculate detailed footprints
#-------------------------------------------------------------------------
extensions <- colnames(E)[7:10]
consumption_categories <- Y.codes$`Final Demand Category`[1:7]
countries <- c("US","CA","AU","EU")

country <- "AT"
country <- "EU27"
# calculate footprints
# for(country in countries){
  for(extension in extensions[-2]){
    data <- footprint(country = country, extension = extension, allocation = allocation)
  }
# }
