##--------------------------------------------------------------------------------------------------
##
## Project: Can trial independent and dependent variables predict experimental outcomes
##
## Hypothesis:  Alterations in locomotion, not in Indi_1, associated with experimental treatments 
##              predicts experimental outcomes.  
## 
##--------------------------------------------------------------------------------------------------

## NOTE1 TO USER - This program requires ~ 35 mins of wall-clock time to execute on a single processor 
## (2.3 GHz, 2 Gb RAM/proc), if all questions below are answered "Y". 

## NOTE2 TO USER - Project's use of unpublished experimental data has forced me to anonymize the data (names) and
## context.
##
##--------------------------------------------------------------------------------------------------


# starting with a blank slate
rm(list=ls())
gc()

#----Load Libraries----
library(logging)
library(data.table)
library(ggplot2)
library(dplyr)
library(glmnet)
library(nnet)
library(caret)


#----Sourcing all the function calls-----
source('/Users/sundars/Documents/Research/DataScience_CertificateProgram/Course2_Methods4DataAnalysis/Project/OpenSource/FN_ProjectSS.R')

#---- starting the main program
if (interactive())
{
  basicConfig()
  addHandler(writeToFile, logger="datalogger", file="ProjectSS.log")
  addHandler(writeToFile, logger="Reg_logger", file="ProjectSS_Regression.log")
  addHandler(writeToFile, logger="Reg_optima", file="ProjectSS_Regression_Optima.log")
  addHandler(writeToFile, logger="Category_logger", file="ProjectSS_Category.log")
  addHandler(writeToFile, logger="Pred_logger", file="ProjectSS_Prediction.log")
  
  # Logging information
  loginfo("Setting wd and loading data.", logger="datalogger")
  # setting working directory - assumes that this R code and the data file exist in the working directory
  setwd('/Users/sundars/Documents/Research/DataScience_CertificateProgram/Course2_Methods4DataAnalysis/Project/OpenSource')
  
  
  # reading/cleaning the data 
  datafile <- 'ProjectData_SSF.csv'
  All_data <- load_clean_data(datafile)
  
  data_dt<-as.data.table(All_data)
  
  #---------------------------------------------------------------------------------------------------------------------------------------
  #---- Exploring Patterns
  DT_notrial<-unique(data_dt[,2:9,with=FALSE]) # these measures are independent of TrialNo 
  PlotKey<-readline("Do you want mean sd data plots for all Trial Independent variates (i.e., Indi_1, Indi_2), i.e., Y or N?  ")
  if (PlotKey == "Y")
  {
    a<-Exploratory_Plotting(DT_notrial,,datalogger)
  }
  
  #   
  PlotKey2<-readline("Do you want mean sd data plots for all Trial Dependent variates (e.g., 1 Hz, 2 Hz components), i.e., Y or N?  ")
  meas<-levels(data_dt$Measure) # these are the different measures representative of locomotion
  if (PlotKey2 == "Y")
  {
    for (i in 1:length(meas))
    {
      Test_measure<-meas[i]
      DT_trial<-unique(data_dt[Measure==Test_measure,c(2:dim(data_dt)[2]),with=FALSE]) # these measures are dependent upon TrialNo
      a<-Exploratory_Plotting(DT_trial,Test_measure,datalogger)
    }
  }  
  
  
  #---------------------------------------------------------------------------------------------------------------------------------------  
  #---- determining the measure (i.e., area, CMx, CMy, length, height) of locomotion that can best classifies treatment by type
  meas_2anal<-CategoricalMeasure(data_dt,Category_logger)
  
  
  #---------------------------------------------------------------------------------------------------------------------------------------  
  #----Linear regressions for exploration of patterns in data
  Notrial_Key<-readline('Do you want to perform exploratory linear regressions for the Trial Independent variates (i.e., Indi_1, Indi_2), i.e., Y or N? ')
  if (Notrial_Key == "Y"){Linear_Regression(DT_notrial,,Reg_logger,Reg_optima)}
  
  
  Trial_Key<-readline('Do you want to perform exploratory linear regressions for the Trial Dependent variates (e.g., f_1 Hz, f_2 Hz), i.e., Y or N? ')
  if (Trial_Key == "Y")
  {
    Test_measure<-meas_2anal
    DT_trial<-unique(data_dt[Measure == Test_measure,c(2:dim(data_dt)[2]),with=FALSE])
    Linear_Regression(DT_trial,Test_measure,Reg_logger,Reg_optima)
  }
  
  
  #---------------------------------------------------------------------------------------------------------------------------------------  
  # ----Linear regressions to determine if experimental outcomes can be predicted
  Pred_Notrial_Key<-readline('Do you want to examine if Trial Independent variates (i.e., Indi_1, Indi_2) can predicted outcomes, i.e., Y or N? ')
  if (Pred_Notrial_Key == "Y")
  {
    # performing a 70-30 split, where 70% is used to determine a predictive model via cross-validation
    # and the remaining data (30%) is used in a final test of model predictions
    n_row<-dim(DT_notrial)[1]
    Train_ind<-sample(n_row,n_row*0.7)
    Train<-DT_notrial[Train_ind,]
    Test<-DT_notrial[!Train_ind,]
    Key<-1
    r2_Notrial<-Predictive_Regression(Train,Test,Key,Pred_logger)
  }
  
  Pred_trial_Key<-readline('Do you want to examine if Trial dependent variates (i.e., locomotion) can predicted outcomes, i.e., Y or N? ')
  if (Pred_trial_Key == "Y")
  {
    # performing a 70-30 split, where 70% is used to determine a predictive model via cross-validation
    # and the remaining data (30%) is used in a final test of model predictions
    Test_measure<-meas_2anal
    DT_trial<-unique(data_dt[Measure == Test_measure,c(2:dim(data_dt)[2]),with=FALSE])
    n_row<-dim(DT_trial)[1]
    Train_ind<-sample(n_row,n_row*0.7)
    Train<-DT_trial[Train_ind,]
    Test<-DT_trial[!Train_ind,]
    Key<-2
    r2_Trial<-Predictive_Regression(Train,Test,Key,Pred_logger)
  }
}


