##-----------------------------------------------------------
##
## Functions called by Main_ProjectSS.R
##
##-----------------------------------------------------------



## Function to read/clean data
load_clean_data <- function(datafile,logger=NA)
{
  data_out = read.csv(datafile,stringsAsFactors = FALSE)
  if (nrow(data_out)==0)
  {
    loginfo('error loading data',logger="datalogger")
  } else {
    loginfo('data loaded good',logger="datalogger")
  }
  
  
  # dropping incomplete cases
  ncomp<-complete.cases(data_out)
  if (sum(!ncomp)>0) # i.e., if there are some incomplete cases
  {data_out<-data[ncomp_out,]}
  
  # changing data types prior to analysis
  data_out<-Change_type(data_out)
  del_indices<-11 # removing DC_Offset term
  data_out<-delete_features(data_out,del_indices)
  return(data_out)
}

## Function to delete unwanted variables/features/columns in the data frame
delete_features = function(Data,del_indices,logger=NA)
{
  orig_length<-ncol(Data)
  
  # deleting select cols if data is a data.table, else, if it is a data frame
  if (is.data.table(Data)){Data_out<-Data[,!del_indices,with=FALSE]} else
  {
    datanames<-names(Data)
    datadrop<-datanames[del_indices]
    Data_out<-Data[,!names(Data) %in%datadrop]
  }
  
  if ((ncol(Data_out)+length(del_indices))==orig_length)
  {
    loginfo('successfully removed unnecessary data', logger="datalogger")
  } else
  {
    loginfo('ERROR removing unnecessary data', logger="datalogger")
  }
  return(Data_out)
}


## Function to change data types appropriately
Change_type <- function(data,logger=NA)
{
  data[,1]<-as.factor(data[,1])
  data[,2]<-as.factor(data[,2])
  data[,3]<-as.factor(data[,3])
  data$Day<-as.numeric(unlist(regmatches(data$Day, gregexpr("[0-9]+", data$Day))))
  loginfo('Changed data types appropriately prior to analysis', logger="datalogger")
  return(data)
}

## Function to determine/observe variation of data over Days and Treatment_Type
Exploratory_Plotting <-function(DT,TM=NA,logger=NA)
{
  if (is.na(TM)) {Key<-0} else {Key<-1}
  #   print(paste('The Key value is:',Key))
  if (Key == 0)
  {
    istart<-3 # index of label in DT for use later
    SubDir1<-"Figures"
    if (dir.exists(SubDir1)){loginfo("Figures directory exists",logger="datalogger")}else{
      dir.create(SubDir1, showWarnings = TRUE)
    }
  } else {
    istart<-9 # index of label in DT for use later
    # repeating this in case Key == 0 was never executed!
    SubDir1<-"Figures"
    if (dir.exists(SubDir1)){loginfo(paste(SubDir1,"directory exists"),logger="datalogger")}else{
      dir.create(SubDir1, showWarnings = TRUE)
      loginfo(paste(SubDir1,"directory was newly created"),logger="datalogger")
    }
    SubDir2<-paste(SubDir1,"/",TM,sep="")
    if (dir.exists(SubDir2)){loginfo(paste(SubDir2,"directory exists"),logger="datalogger")}else{
      dir.create(SubDir2, showWarnings = TRUE)
      loginfo(paste(SubDir2,"directory was newly created"),logger="datalogger")
    }
  }
  
  if (Key == 0){del_indices<-3} else {del_indices<-c(3:9)}
  
  mean_dt<-DT[,lapply(.SD,mean),by=list(Day,Treatment_Type)]
  mean_dt<-delete_features(mean_dt,del_indices)
  sd_dt<-DT[,lapply(.SD,sd),by=list(Day,Treatment_Type)]
  sd_dt<-delete_features(sd_dt,del_indices)
  
  MeSD<-left_join(mean_dt,sd_dt, by=c("Treatment_Type","Day"))
  nvar<-length(names(MeSD))-2
  
  sd_max<-MeSD[,3:(nvar/2+2),with=FALSE]+MeSD[,(nvar/2+3):(nvar+2),with=FALSE]
  sd_min<-MeSD[,3:(nvar/2+2),with=FALSE]-MeSD[,(nvar/2+3):(nvar+2),with=FALSE]
  
  for (i in 1:(nvar/2))
  {
    
    DF<-data.frame(MeSD[,2,with=FALSE],MeSD[,1,with=FALSE],MeSD[,i+2,with=FALSE],sd_min[,i,with=FALSE],sd_max[,i,with=FALSE])
    names(DF)<-c('Treattype','x1','y1','sd1','sd2')
    
    ylab<-names(DT)[i+istart]
    xlab<-names(DT)[3]
    if (Key == 0)
    {
      title1<-paste(ylab,"change over time")
      file1<-paste(SubDir1,"/",ylab,"_fn_Treatment.jpg",sep="")
    } else 
    {
      title1<-paste(TM,":", ylab,"change over time")
      file1<-paste(SubDir2,"/",ylab,"_fn_Treatment.jpg",sep="")
    }
    
    graphics.off()
    
    gg<-ggplot(DF,aes(x=x1,y=y1))+geom_point()+geom_line()+facet_wrap(~Treattype,ncol=5)+
      geom_errorbar(aes(ymin=sd1, ymax=sd2))+
      labs(x=xlab,y=ylab,title=title1) +theme_bw()
    
    ggsave(filename=file1,plot=gg)
    
    loginfo(paste("plotted",file1,"successfully"), logger="datalogger")
  }
  return(NULL)
}

## Function to determine the gait trial that is most representative of mean gait features for a given Treatment type
# Match_Features identifies the DT index at which the fft components of individual gait trials are closest to the 
# mean of the trial based on eucledian distance on a n-dimensional space (where n =  no of fft components x no. of days)
Match_Features_Index<-function(DT_trial)
{
  mean_dt<-DT_trial[,lapply(.SD,mean),by=list(Spec_ID,Treatment_Type,Day)]
  key1<-names(DT_trial)
  setkeyv(DT_trial,cols=key1)
  setkeyv(mean_dt,cols=key1)
  DT_trial$Index<-1:dim(DT_trial)[1]# setting the index so the Match_index can be computed
  
  nmean<-dim(mean_dt)[1]
  Match_index<-rep(NA,nmean)
  
  for (i in 1:nmean)
  {
    junk_DT<-DT_trial[((Spec_ID == mean_dt[i,Spec_ID]) & (Treatment_Type == mean_dt[i,Treatment_Type]) & (Day == mean_dt[i,Day])), ]
    nt<-dim(junk_DT)# no of trials in this particular cohort
    junk2_DT<-rbind(mean_dt[i,],junk_DT[,1:(nt[2]-1),with=FALSE])
    junk2_DT<-junk2_DT[,c(5:dim(junk2_DT)[2]),with=FALSE]
    dist2_pt<-dist(junk2_DT)
    Match_index[i]<-junk_DT[which.min(dist2_pt[1:nt[1]]),Index] # trial index closest to the mean gait features
  }
  return(Match_index)
}

## Function to perform joining operation across days, and appropriately truncating/naming the data table
DayJoin<-function(DT)
{
  uniq_days<-unique(DT$Day)
  ndays<-length(uniq_days)
  DT_join<-DT[Day == uniq_days[1],]
  DT_join[,Day:=NULL]
  #   DT_join[,TrialNo:=NULL]
  new_names<-sapply(5:dim(DT)[2],function(k){paste(names(DT)[k],"_d",uniq_days[1],sep="")})
  if (ndays > 1)
  {
    for (j in 2:ndays)
    {
      DT_join<-merge(DT_join,DT[Day == uniq_days[j]],by=c("Spec_ID","Treatment_Type","TrialNo"))
      DT_join[,Day:=NULL]
      new_names<-c(new_names,sapply(5:dim(DT)[2],function(k){paste(names(DT)[k],"_d",uniq_days[j],sep="")}))
    }
  }
  new_names<-c("Spec_ID","Treatment_Type",new_names)
  DT_join[,TrialNo:=NULL]
  names(DT_join)<-new_names
  DT_join[,Spec_ID:=NULL]
  return(DT_join)
}

## Function to determine which gait measure most categorically represents Treatment types.
# Classification by multinomial logistic regression
CategoricalMeasure<- function(DT_Orig,logger=Category_logger)
{
  meas<-levels(DT_Orig$Measure)
  best_meas<-meas[1]
  best_acc<-0
  
  for (i in 1:length(meas))
  {
    DT_M<-DT_Orig[Measure == meas[i]]
    del_indices<-c(1,5:9)
    DT<-delete_features(DT_M,del_indices)
    DT<-data.table(DT[,c(1:4),with=FALSE],scale(DT[,c(5:dim(DT)[2]),with=FALSE]))# scaling the gait fft components
    
    Feat_Index<-Match_Features_Index(DT)# extracts the trial that 'best' represent the mean of the fft components of gait
    DT_test<-DT[!Feat_Index,]
    DT<-DT[Feat_Index,]
    
    DT$TrialNo<-1 #assigning a dummy trial no for the extracted trial
    DT_join<-DayJoin(DT)
    
    predictor<-DT_join[,Treatment_Type]
    xmod<-model.matrix(Treatment_Type ~ .,DT_join)[,-1]
    dt_lasso<-glmnet(xmod,predictor,alpha=1,family='multinomial',type.multinomial="grouped") # multinomial lasso regression
    plot(dt_lasso,xvar='lambda',label=TRUE)
    
    dt_cv<-cv.glmnet(xmod,predictor,alpha=1,family='multinomial',type.multinomial="grouped") # cross-validation
    plot(dt_cv)
    
    lambda_min<-dt_cv$lambda.min
    best_coef<-lapply(coef(dt_lasso),function(x) x[,dt_lasso$lambda == lambda_min])# determining the best coeffs for the models for each Treatment type factor
    coef_gt0<-lapply(best_coef,function(x) which(abs(x)>1e-8))# determining the id of coeffs that are not negligible
    
    uniq_gt0_coef<-unique(unlist(coef_gt0, use.names = FALSE)) # finding unique non-zero coeffs for each of the models
    
    #now, that we have a multinomial logistic regression model form (i.e., independent variates) identified as potentially predictive
    # examining its predictive utility
    DT_join$Treattype<-relevel(DT_join$Treatment_Type,ref="Treatment1")
    formula<-as.formula(paste("Treattype ~",paste0(names(DT_join)[uniq_gt0_coef[-1]],collapse="+",sep=""))) # formula for 'optimal' model
    mult_logreg<-multinom(formula,DT_join)# multinomial logistic regression using 'optimal' model    
    fit_onTrain<-fitted(mult_logreg)# making sure that the model can predict the data it was trained on.
    predict_Train<-apply(fit_onTrain,1,which.max)
    CM_train<-confusionMatrix(predict_Train,as.numeric(DT_join$Treattype))
    loginfo(paste("The accuracy of the measure:",meas[i], ", in classifying Treatment_Type from the TRAINING data is:",CM_train$overall[1]),logger="Category_logger")
    
    # obtaining the test data for the model, and predicting the Test data
    DT_test_join<-DayJoin(DT_test)
    DT_test_join$Treattype<-relevel(DT_test_join$Treatment_Type,ref="Treatment1")
    predict_Test<-apply(predict(mult_logreg,DT_test_join,"probs"),1,which.max)
    CM_Test<-confusionMatrix(predict_Test,as.numeric(DT_test_join$Treattype))# confusion matrix for multinomial classifier
    mod_acc<-signif(CM_Test$overall[1],3)
    loginfo(paste("The accuracy of the measure:",meas[i], ", in classifying Treatment_Type from the TEST data is:",mod_acc),logger="Category_logger")
    if (best_acc < mod_acc)
    {
      best_acc<-mod_acc
      best_meas<-meas[i]
    }
  }
  loginfo(paste("The gait measure that best classifies Treatment type is:",best_meas))
  return(best_meas)
}

## Function to determine p value from a linear regression model 'fit'
lm_p <- function(fit) {
  if (class(fit) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(fit)$fstatistic
  if (length(f) == 0) # i.e., the f-statistic is not available
  {
    p <- 1
  } else {p <- pf(f[1],f[2],f[3],lower.tail=F)}
  
  attributes(p) <- NULL
  return(p)
}

## Function to perform linear regressions for the purpose of exploration
Linear_Regression<-function(DT,TM=NA,logger=Reg_logger,logger2=Reg_optima)
{
  if (is.na(TM)) # then working on notrial data
  {tit<-"on trial independent data"
   Key<-0} else{
     tit<-paste("on trial dependent",TM,"data")
     Key<-1
   }
  
  # this code is to make sure that the file directories exist so that plots can be saved appropriately
  dirA<-"Figures_Regression"
  if (dir.exists(dirA)){loginfo(paste(dirA,"directory exists"),logger="Reg_logger")} else{
    dir.create(dirA, showWarnings = TRUE)
    loginfo(paste(dirA,"directory was newly created"),logger="Reg_logger")
  }
  if (Key == 0)
  {
    SubDir1<-paste(dirA,'/No_trial_data',sep="")
    SubDir2<-SubDir1
    if (dir.exists(SubDir1)){loginfo(paste(SubDir1,"directory exists"),logger="Reg_logger")}else{
      dir.create(SubDir1, showWarnings = TRUE)
      loginfo(paste(SubDir1,"directory was newly created"),logger="Reg_logger")
    }
  } else {
    # repeating this in case Key == 0 was never executed!
    SubDir1<-"Figures_Regression/Trial_data"
    if (dir.exists(SubDir1)){loginfo(paste(SubDir1,"directory exists"),logger="Reg_logger")}else{
      dir.create(SubDir1, showWarnings = TRUE)
      loginfo(paste(SubDir1,"directory was newly created"),logger="Reg_logger")
    }
    SubDir2<-paste(SubDir1,"/",TM,sep="")
    if (dir.exists(SubDir2)){loginfo(paste(SubDir2,"directory exists"),logger="Reg_logger")}else{
      dir.create(SubDir2, showWarnings = TRUE)
      loginfo(paste(SubDir2,"directory was newly created"),logger="Reg_logger")
    }
  }
  
  loginfo("######################################################################",logger="Reg_logger")
  loginfo(paste("Performing linear regressions", tit),logger="Reg_logger")
  loginfo("######################################################################",logger="Reg_optima")
  loginfo(paste("Performing linear regressions", tit),logger="Reg_optima")
  
  n_dep<-3 # i.e., dependent variates are Outcome_1, Outcome_2, Outcome_3
  if (Key == 0)# i.e. for trial independent data
  {
    n_indi<-2 #i.e., the independent variables are: Indi_1, Indi_2
    istart<-3
  } else {
    n_indi<-30 #i.e., independent variables are all the fft components of gait
    istart<-9
  }
  for (j in 1:n_dep)
  {
    loginfo("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^",logger="Reg_logger")
    loginfo("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^",logger="Reg_optima")
    for (i in 1:n_indi)
    {
      loginfo("------------------------------------------------------------------------",logger="Reg_logger")
      loginfo("------------------------------------------------------------------------",logger="Reg_optima")
      
      if (Key == 1){loginfo(paste("For the measure:",TM),logger="Reg_logger")}
      loginfo(paste("Linear regression of",names(DT)[5+j],"vs",names(DT)[istart+i]),logger="Reg_logger")
      
      # four models: single factor, single factor & Treatment type, single factor & day, and single factor, day & Treatment type
      for (k in 1:4)
      {
        if (k == 1)# single factor
        {
          formula<-as.formula(paste(names(DT)[5+j], "~", names(DT)[istart+i]))
          loginfo(paste("Regression on:",names(DT)[istart+i],logger="Reg_logger"))
          SubDir3<-paste(SubDir2,"/SingleFactor",sep="")
          if (dir.exists(SubDir3)){loginfo(paste(SubDir3,"directory exists"),logger="Reg_logger")}else
          {
            dir.create(SubDir3, showWarnings = TRUE)
            loginfo(paste(SubDir3,"directory was newly created"),logger="Reg_logger")
          }
        }
        if (k == 2) # single factor & Treatment type
        {
          formula<-as.formula(paste(names(DT)[5+j],"~",names(DT)[istart+i],"+",names(DT)[2],"+",names(DT)[istart+i],"*",names(DT)[2]))
          loginfo(paste("Regression on Treatment type and:",names(DT)[istart+i]),logger="Reg_logger")
          SubDir3<-paste(SubDir2,"/SingleFactor_Treatment_Type",sep="")
          if (dir.exists(SubDir3)){loginfo(paste(SubDir3,"directory exists"),logger="Reg_logger")}else
          {
            dir.create(SubDir3, showWarnings = TRUE)
            loginfo(paste(SubDir3,"directory was newly created"),logger="Reg_logger")
          }                  
        }
        if (k == 3)# single factors & Day
        {
          formula<-as.formula(paste(names(DT)[5+j],"~",names(DT)[istart+i],"+",names(DT)[3],"+",names(DT)[istart+i],"*",names(DT)[3]))
          loginfo(paste("Regression on Day and:",names(DT)[istart+i]),logger="Reg_logger")
          SubDir3<-paste(SubDir2,"/SingleFactor_Day",sep="")
          if (dir.exists(SubDir3)){loginfo(paste(SubDir3,"directory exists"),logger="Reg_logger")}else
          {
            dir.create(SubDir3, showWarnings = TRUE)
            loginfo(paste(SubDir3,"directory was newly created"),logger="Reg_logger")
          }
        }
        if (k == 4) # single factor, Treatment type and day
        {
          formula<-as.formula(paste(names(DT)[5+j],"~",names(DT)[istart+i],"+",names(DT)[2],"+",names(DT)[3],
                                    "+",names(DT)[istart+i],"*",names(DT)[2],"+",names(DT)[istart+i],"*",names(DT)[3],
                                    "+",names(DT)[2],"*",names(DT)[3],
                                    "+",names(DT)[istart+i],"*",names(DT)[2],"*",names(DT)[3]))
          loginfo(paste("Regression on Treatment type, day and:",names(DT)[istart+i]),logger="Reg_logger")
          SubDir3<-paste(SubDir2,"/SingleFactor_Treatment_Type_Day",sep="")
          if (dir.exists(SubDir3)){loginfo(paste(SubDir3,"directory exists"),logger="Reg_logger")}else
          {
            dir.create(SubDir3, showWarnings = TRUE)
            loginfo(paste(SubDir3,"directory was newly created"),logger="Reg_logger")
          }
        }
        
        fit_pre<-lm(formula,DT) # linear regression
        if (Key == 0){fit<-step(fit_pre)} else {fit<-fit_pre}# determining optimal model per AIC for no trial data (optima for trial data determined later)
        pval<-signif(lm_p(fit),3)
        r2<-signif(summary(fit)$r.squared,3)
        r2a<-signif(summary(fit)$adj.r.squared,3)
        loginfo(paste("the R square =",r2," adjusted R square =",r2a," and the p value =",pval),logger="Reg_logger")
        if (Key == 0)
        {
          loginfo(paste("Linear regression of",names(DT)[5+j],"vs",names(DT)[istart+i]),logger="Reg_optima")
          loginfo("The optimal model, per AIC, for the no-trial data has the following form:",logger="Reg_optima")
          loginfo(fit$call$formula,logger="Reg_optima")
          loginfo(paste("the R square =",r2," adjusted R square =",r2a," and the p value =",pval),logger="Reg_optima")
        }
        fit_res<-lm(fit$fitted.values ~ fit$residuals)
        pval_res<-signif(lm_p(fit_res),3)
        if (pval_res <=0.05){logwarn('The error is not normally distributed, and the variates may require transformation', logger="Reg_logger")}else{
          loginfo('The error is normally distributed',logger="Reg_logger")
        }
        
        # plotting the data and the regression line for trial independent variates and
        # for trail dependent variates when p <= 0.05
        if ((Key == 0) | ((Key == 1) & (pval <= 0.05)))
        {
          graphics.off()
          xlab<-names(DT)[istart+i]
          ylab<-names(DT)[5+j]
          file1<-paste(SubDir3,"/",ylab,"_(vs)_",xlab,".jpg",sep="")
          if (Key == 1){file1<-paste(SubDir3,"/",ylab,"_(vs)_",xlab,".jpg",sep="")}
          
          DF<-data.frame(DT[,2,with=FALSE],DT[,3,with=FALSE],DT[,i+istart,with=FALSE],DT[,j+5,with=FALSE])
          names(DF)<-c("Treattype",'Day','x','y')
          
          if (k == 1) # single factor
          {
            title1<-paste("Single Factor Regression","\n",names(DT)[5+j],"vs",names(DT)[istart+i],"\nadj r square =",r2a," and  p =",pval)
            gg<-ggplot(DF,aes(x,y))+geom_point()+geom_smooth(method=lm,fullrange=TRUE)+
              labs(x=xlab,y=ylab,title=title1)+theme_bw()
          }
          if (k == 2) # single factor and Treatment type
          {
            title1<-paste("Treatment_Type & Individual Gait Component Regression","\n",names(DT)[5+j],"vs",names(DT)[istart+i],"\nadj r square =",r2a," and p =",pval)
            gg<-ggplot(DF,aes(x,y))+geom_point()+geom_smooth(method=lm,fullrange=TRUE)+
              labs(x=xlab,y=ylab,title=title1)+theme_bw()+facet_wrap(~Treattype)
          }
          if (k == 3) # single factor and day
          {
            title1<-paste("Day & Individual Gait Component Regression","\n",names(DT)[5+j],"vs",names(DT)[istart+i],"\nadj r square =",r2a," and p =",pval)
            gg<-ggplot(DF,aes(x,y))+geom_point()+geom_smooth(method=lm,fullrange=TRUE)+
              labs(x=xlab,y=ylab,title=title1)+theme_bw()+facet_wrap(~Day)
          }
          if (k == 4) # single factor, Treatment type and day
          {
            title1<-paste("Treatment_Type, Day and Individual Component Regression","\n",names(DT)[5+j],"vs",names(DT)[istart+i],"\nadj r square =",r2a," and p =",pval)
            DF$Day<-as.factor(DF$Day)
            gg<-ggplot(DF,aes(x,y))+geom_point(aes(color=Day))+geom_smooth(method=lm,fullrange=TRUE)+
              labs(x=xlab,y=ylab,title=title1)+theme_bw()+facet_wrap(~Treattype)
          }
          ggsave(filename=file1,plot=gg)
        }
      } # k index of for statement; type of regression: single factor or 1st, 2nd order
    } # i index of for statement; independent vars
    
    # regression of all gait components combined for the trial dependent data, and det. of optimal model per AIC
    if (Key == 1)
    {
      # doing single factor, single factor & Treatment type, single factor & day, and single factor, day & Treatment type
      loginfo(paste("Linear regression of:",names(DT)[5+j],"vs all gait components"),logger="Reg_optima")
      
      for (k in 1:5)
      {
        if (k == 1) # all gait components
        {
          DT2<-DT[,c(5+j,10:dim(DT)[2]),with=FALSE]
          formula<-as.formula(paste(names(DT2)[1],"~ ."))
          loginfo("Analysis of Gait components",logger="Reg_optima")
        }
        if (k == 2) # all gait components and Treatment type
        {
          DT2<-DT[,c(2,5+j,10:dim(DT)[2]),with=FALSE]
          formula<-as.formula(paste(names(DT2)[2],"~ . +(. - Treatment_Type)*Treatment_Type"))
          loginfo("Analysis of Gait components and Treatment_Type",logger="Reg_optima")
        }
        if (k == 3) # all gait components and day
        {
          DT2<-DT[,c(3,5+j,10:dim(DT)[2]),with=FALSE]
          formula<-as.formula(paste(names(DT2)[2],"~ . +(. - Day)*Day"))
          loginfo("Analysis of Gait components and Day",logger="Reg_optima")
          
        }
        if (k == 4) # all gait components, Treatment type and day
        {
          DT2<-DT[,c(2,3,5+j,10:dim(DT)[2]),with=FALSE]
          formula<-as.formula(paste(names(DT2)[3],"~ . + Treatment_Type*Day + (. - Treatment_Type - Day)*Treatment_Type",
                                    "+ (. - Treatment_Type - Day)*Day + (. - Treatment_Type - Day)*Treatment_Type*Day"))
          loginfo("Analysis of Gait components, Treatment_Type, and Day",logger="Reg_optima")
          
        }
        if (k == 5) # all gait components, Indi_1, Indi_2, Treatment type and day <- the FULL regression
        {
          DT2<-DT[,c(2:5,5+j,10:dim(DT)[2]),with=FALSE]
          formula<-as.formula(paste(names(DT2)[5],"~ . + Treatment_Type*Day + (. - Treatment_Type - Day)*Treatment_Type",
                                    "+ (. - Treatment_Type - Day)*Day + (. - Treatment_Type - Day)*Treatment_Type*Day"))
        }
        fit_pre<-lm(formula,DT2) # linear regression
        if (is.nan(summary(fit_pre)$adj.r.squared))
        {
          loginfo(paste("the current model is bizzare and provides a perfect fit, with no residual errors for",names(DT[5+j])),logger="Reg_optima")
        } else
        {
          fit<-step(fit_pre) # optimal model per min AIC
          pval<-signif(lm_p(fit),3)
          r2<-signif(summary(fit)$r.squared,3)
          r2a<-signif(summary(fit)$adj.r.squared,3)
          loginfo("The optimal model takes the following form",logger="Reg_optima")
          loginfo(fit$call$formula,logger="Reg_optima")
          loginfo(paste("the R square =",r2,"the adjusted R square=",r2a,"and the p -value = ",pval),logger="Reg_optima")
        }
      }
    }
  } # j index of for statement; dependent vars
  return(NULL)
}

## Function to perform linear regression for the purpose of prediction
Predictive_Regression<-function(DT_Train,DT_Test,Key,Pred_logger)
{
  # this code is to make sure that the file directories exist so that plots can be saved appropriately
  dirA<-"Figures_Prediction"
  if (dir.exists(dirA)){loginfo(paste(dirA,"directory exists"),logger="Pred_logger")} else{
    dir.create(dirA, showWarnings = TRUE)
    loginfo(paste(dirA,"directory was newly created"),logger="Pred_logger")
  }

  # for trial independent data set (i.e., Key=1), three models (i.e., Indi_1, Indi_2, Indi_1+Indi_2)xDay
  # for trial dependent data set (i.e., Key=2), four models (i.e., Gait, GaitxIndi_1, GaitxIndi_2, GaitxIndi_1xIndi_2)xDay
  # in these cases, Treatment_Type is ignored as a factor, because given gait, Indi_1 and Indi_2 across time, the idea
  # is to develop a model that can ultimately predict Outcomes for a new Treatment type that has not yet been implemented
  DT_Train<-DT_Train[,c(3:dim(DT_Train)[2]),with=FALSE]
  DT_Test<-DT_Test[,c(3:dim(DT_Test)[2]),with=FALSE]
  
  n_dep<-3 # dependent variates are Outcome_1, Outcome_2, and Outcome_3
  
  loginfo("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$",logger="Pred_logger")
  if (Key == 1){TnoT<-"No_Trial"}else{TnoT<-"Trial"}
  loginfo(paste("Performing predictive analysis on",TnoT,"Data"),logger="Pred_logger")
  
  n_mods<-5
  
  r2_mat<-matrix(NA,nrow = n_dep,ncol=n_mods*2)
  
  for (i in 1:n_dep) 
  {
    best_r2<-0
    best_mod<-0
    loginfo("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^",logger="Pred_logger")
    loginfo(paste("For the dependent variable: ",names(DT_Train)[3+i]),logger="Pred_logger")
    
    for (j in 1:n_mods)
    {
      loginfo("---------------------------------------------------------------------",logger="Pred_logger")
      loginfo(paste("For Model #",j),logger="Pred_logger")
      if (Key == 1)# no trial data, creating 3 models
      {
        if (j == 1) # Indi_1
        {
          mod_train<-DT_Train[,c(2,3+i),with=FALSE]
        }
        if (j == 2) # Indi_2
        {
          mod_train<-DT_Train[,c(3,3+i),with=FALSE]
        }
        if (j == 3) # Indi_1, Day
        {
          mod_train<-DT_Train[,c(1,2,3+i),with=FALSE]
        }
        if (j == 4)# Indi_2, Day
        {
          mod_train<-DT_Train[,c(1,3,3+i),with=FALSE]
        }
        if (j == 5)# Indi_1, Indi_2, Day
        {
          mod_train<-DT_Train[,c(1:3,3+i),with=FALSE]
        }        
      }
      if (Key == 2)# Trial dependent data, creating 4 models
      {
        if (j == 1) # Gait
        {
          mod_train<-DT_Train[,c(8:dim(DT_Train)[2],3+i),with=FALSE]
        }
        if (j == 2)# Gait, Day
        {
          mod_train<-DT_Train[,c(1,8:dim(DT_Train)[2],3+i),with=FALSE]
        }
        if (j == 3)# Gait,Indi_1,Day
        {
          mod_train<-DT_Train[,c(1,2,8:dim(DT_Train)[2],3+i),with=FALSE]
        }
        if (j == 4)# Gait, Indi_2, Day
        {
          mod_train<-DT_Train[,c(1,3,8:dim(DT_Train)[2],3+i),with=FALSE]
        }
        if (j == 5) # Gain, Indi_1, Indi_2, Day
        {
          mod_train<-DT_Train[,c(1:3,8:dim(DT_Train)[2],3+i),with=FALSE]
        }
      }
      if ((Key == 1 & j < 3) | (Key == 2 & j == 1))
      {
        formula<-as.formula(paste(names(DT_Train)[3+i],"~ ."))
      } else {formula<-as.formula(paste(names(DT_Train)[3+i],"~ . + Day + (. - Day)*Day"))}
      
      fullMod_fit<-lm(formula,mod_train)# All terms model
      pred_fullmod<-predict(fullMod_fit,DT_Test)
      predVSact<-data.frame(pred_fullmod,DT_Test[,3+i,with=FALSE])
      names(predVSact)<-c("Predicted","Actual")
      fit2_fullmod<-lm(Actual ~ Predicted, predVSact)
      r2<-signif(summary(fit2_fullmod)$r.squared,3)
      r2a<-signif(summary(fit2_fullmod)$adj.r.squared,3)
      pval<-signif(lm_p(fit2_fullmod),3)
      r2_mat[i,j]<-r2
      title1<-paste(names(DT_Train)[3+i],"~",names(mod_train)[1])
      title2<-paste(title1,"\n","adj r square =",r2a, "and p value = ",pval)
      xlab<-paste("Predicted",names(DT_Train)[3+i])
      ylab<-paste("Actual",names(DT_Train)[3+i])
      loginfo(paste("For the full model, the rsquare =",r2,"the adj r square = ",r2a,"and the p-value =",pval),logger="Pred_logger")
      
      # determining an 'optimal' model via lasso regression and cross-validation
      if ((Key == 1 & j > 2) | (Key == 2 & j > 1)) # only for models 1, 2 for No_trial data and model 1 for trial data
      {
        predictor<-unlist(mod_train[,dim(mod_train)[2],with=FALSE])
        xmod<-model.matrix(formula,mod_train)[,-1]
        dt_lasso<-glmnet(xmod,predictor,alpha=1,family='gaussian')
        plot(dt_lasso,xvar='lambda',label=TRUE)
        dt_cv<-cv.glmnet(xmod,predictor,alpha=1,family='gaussian')
        plot(dt_cv)
        lambda_min<-dt_cv$lambda.min
        best_coef<-coef(dt_lasso)[,dt_lasso$lambda == lambda_min]
        coef_gt0<-which(abs(best_coef)>1.e-8)# determining the id of coeffs that are not negligible
        
        #now, that we have a regression model form (i.e., independent variates) identified as potentially predictive
        # examining its predictive utility
        xmod_dt<-data.table(xmod)
        title1<-paste(names(DT_Train)[3+i]," ~ ",paste0(names(xmod_dt)[coef_gt0[-1]-1],collapse="+",sep=""))
        loginfo(paste("The best form of the model is:",title1),logger="Pred_logger")
        
        formula_new<-as.formula(title1)
        fit1<-lm(formula_new,DT_Train)
        pred_test<-predict(fit1,DT_Test)
        plot(pred_test,unlist(DT_Test[,3+i,with=FALSE]))
        predVSact<-data.frame(pred_test,DT_Test[,3+i,with=FALSE])
        names(predVSact)<-c("Predicted","Actual")
        fit2<-lm(Actual ~ Predicted, predVSact) # for how good the predictions are
        xlab<-paste("Predicted",names(DT_Train)[3+i])
        ylab<-paste("Actual",names(DT_Train)[3+i])
        r2<-signif(summary(fit2)$r.squared,3)
        r2a<-signif(summary(fit2)$adj.r.squared,3)
        pval<-signif(lm_p(fit2),3)
        r2_mat[i,n_mods+j]<-r2a
        title2<-paste(title1,"\n","adj r square =",r2a, "and p value = ",pval)
        
        loginfo(paste("For this form of model, the r square =",r2, "the adj r square = ",r2a," and the p value=",pval),logger="Pred_logger")
      }
      
      
      graphics.off()
      file1<-paste(dirA,"/",TnoT,"_",names(DT_Train)[3+i],"_Model",j,".jpg",sep="")
      gg<-ggplot(predVSact,aes(Predicted,Actual))+geom_point()+geom_smooth(method=lm,fullrange=TRUE)+
        labs(x=xlab,y=ylab,title=title2)+theme_bw()
      ggsave(filename=file1,plot=gg)
      if (best_r2<r2)
      {
        best_r2<-r2
        best_mod<-j
      }
    }# for (j in 1:n_mod)
    loginfo("---------------------------------------------------------------------",logger="Pred_logger")
    loginfo(paste("For",names(DT_Train)[3+i]," the best predictor is model #",best_mod,"with an r square =",best_r2),logger="Pred_logger")
  }# for (i in 1:n_dep)
  
  return(r2_mat)
}


