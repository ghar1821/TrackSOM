library(Spectre)
library(gtools)
library(Biobase)

setwd("M:/givanna/COVID_timeseries/immune_panel")
data.dir <- getwd()

list.files(data.dir, "*.csv")

data.list <- Spectre::read.files()
cell.dat <- data.table::rbindlist(data.list, fill = TRUE)

## read in metadata
meta.dat <- data.table(read.csv("M:/givanna/COVID_timeseries/metadata/CVD-maser-CSV.csv"))
meta.dat <- meta.dat[Inclusion == 'Include', ]
meta.unique.sample <- unique(meta.dat$SampleIndex)
cell.dat.unique.sample <- unique(cell.dat$SampleIndex)
setdiff(meta.unique.sample, cell.dat.unique.sample)
### Some samples are in meta data but not in dataset. 

## Get only COVID patient and remove those treated at home or NP
cell.dat <- cell.dat[Diagnosis == 'COVID' & DiseaseSeverity %in% c('COVID_ward', 'COVID_ICU', 'COVID_hospital'),]

## Get data for clustering cols
clustering.cols <- names(cell.dat)[61:71]

save.dir <- "data.binned.4tp"
setwd(save.dir)

unique.dpo <- mixedsort(unique(cell.dat$DPO))
## Because we have 41 time point, we just merge the last timepoint into same last bin.
for (x in c(0:9)) {
  idx_start <- (x * 4) + 1
  
  ## Because we have 41 time point, we just merge the last timepoint into same last bin.
  if (x == 9) {
    idx_end <- length(unique.dpo)
  } else {
    idx_end <- (x+1) * 4
  }
  
  filtered.dat.list <- lapply(c(idx_start:idx_end), function(idx) {
    cell.dat[DPO == unique.dpo[idx], ]
  })
  filtered.dat <- rbindlist(filtered.dat.list, fill=TRUE)
  print(unique(filtered.dat$DPO))
  Spectre::write.files(filtered.dat, paste0("AllFeatures_TPBin_", x),
                       write.csv = TRUE,
                       write.fcs = TRUE)
  
  Spectre::write.files(filtered.dat[, ..clustering.cols], paste0("ClustFeatures_TPBin_", x),
                       write.csv = TRUE,
                       write.fcs = TRUE)
}
