##    Rewrite attempt December 13 2018 by Zach Quinlan
# I am trying to make the matlab scirpt I used for fDOM in R. 


# Clear current environemnt
rm(list =ls())
#Clear plots
if(is.null(dev.list())) dev.off()


EEMfilename = 'CSV.CSV'
# Note: column 2 should list text Blank or Sample, column 1 should have unique informative names of samples,
#       columns 3, 4, and 5 should list complete path names for each sample and respective blanks,
#       and column 6 lists complete path names for the absorbance scan.


slit_width <- 5
ex_start <- 240
ex_end <- 500
ex_step <- 5
em_start <-248  # these values are designed to estimate the range of emissions matching the number of Aqualog export rows...
em_end <- 824.6 #...assuming there is a EmStep distance between each row or each emission value.
em_step <- 4.65 #On the aqualog this is dictated by the number of pixels assigned to integrate 4.65nm = 8 pixels

# these develop axis values for Ex/Em, example 500 495 490 ... 245 240
ex <- seq(ex_start, ex_end, by = ex_step)
em <- seq(em_start, em_end, by = em_step)

