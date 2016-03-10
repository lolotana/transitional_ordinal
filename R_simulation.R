
############################################################################################################################
### Autoregressive transitional ordinal model to test for treatment effect in neurological trials with complex endpoints ###
############################################################################################################################

# README:

# 1) This is the simplified & commented R code used for the simulation part of the manuscript 
# 2) Refer to the manuscript for deails on the simulation study
# 3) The R code only runs if provided with patient data formatted accordingly to the given instructions
# 4) The patient data is not publicly available and not shipped with this script, you may apply for data access at emsci.org
# 5) The R code is set up to run on a simplified version of the whole simulation. 
#   Change settings (original settings commented out) under the "set simulation settings" code chunck 
#    to run full simulation. The latter is very (!) time-consuming, you may want to run 
#    it on a multi-core server and/or parallelize it.


### import patient data ##########################################################################

# store your data in an R object called "mydata" !

# Each row represents one of the key muscles in the Upper Extremity Motor Scores, 
# being at or below the Motor Level of the left and right body side, respectively


# a fictitious patient would look like this (first row are column names):

# id    side      Motor_2w    diff    Md        MS_2w   MS_6m   difnum
# 0001  left      C5          0       C5.0      4       5       0
# 0001  left      C5          -1      C5.-1     4       5       -1
# 0001  left      C5          -2      C5.-2     3       5       -2
# 0001  left      C5          -3      C5.-3     2       5       -3
# 0001  left      C5          -4      C5.-4     1       2       -4
# 0001  right     C6          0       C6.0      4       5       0
# 0001  right     C6          -1      C6.-1     3       5       -1
# 0001  right     C6          -2      C6.-2     2       5       -2
# 0001  right     C6          -3      C6.-3     1       3       -3


# colums: 
#     id - patient ID
#     side - body side of the key muscle [left, right]
#     Motor_2w - motor level at 2 weeks/baseline
#     diff - distance of key muscle from  motor level (measured as mumber of key muscle)
#     Md - combination of Motor_2w and diff
#     MS_2w - motor score at 2 weeks/baseline
#     MS_6m - motor score at 6 months/end
#     difnum - numeric version of diff for later manipulations



### R session ################################################################################

set.seed(2709)

library("coin")
library("MASS")


### set simulation settings ###################################################################

# set the number of iterations
nsim <- 2
# nsim <- 1000

# naming of iterations
simid <- factor(paste("sim", formatC(1:nsim, width = nchar(as.character(nsim)), flag = "0"), sep = "_"))

# set total trial size and treatment effect 
totaltrialsize <- seq(100, 200, by=100)
simulatedtreatment <- round(log((10:11)/10), digits = 4)
# totaltrialsize <- seq(50, 200, by=25)
# simulatedtreatment <- round(log((10:15)/10), digits = 4)

# produce combinations and save
param <- expand.grid(simid = simid, n = totaltrialsize, treat =  simulatedtreatment )
rownames(param) <- do.call("paste", c(param, sep = ":"))


### compute autoregressive part and fit reference model ####################################

ref_mod<-function(){

  # select blocks of data with different distance from the lesion level
  level0<-subset(mydata, difnum == 0)
  level1<-subset(mydata, difnum == -1)
  level2<-subset(mydata, difnum == -2)
  level3<-subset(mydata, difnum == -3)
  level4<-subset(mydata, difnum == -4)
  
  # sort blocks by patient ID and side
  a<-subset( merge(level4,level3, by=c("id", "side"), all=FALSE), select= c(id, side, MS_6m.y) )
  aa<- a[order(a$id, a$side),]
  level4$MS_6m_auto<-aa$MS_6m.y
  
  b<-subset( merge(level3,level2, by=c("id", "side"), all=FALSE), select= c(id, side, MS_6m.y) )
  bb<- b[order(b$id, b$side),]
  level3$MS_6m_auto<-bb$MS_6m.y
  
  c<-subset( merge(level2,level1, by=c("id", "side"), all=FALSE), select= c(id, side, MS_6m.y) )
  cc<- c[order(c$id, c$side),]
  level2$MS_6m_auto<-cc$MS_6m.y
  
  d<-subset( merge(level1,level0, by=c("id", "side"), all=FALSE), select= c(id, side, MS_6m.y) )
  dd<- d[order(d$id, d$side),]
  level1$MS_6m_auto<-dd$MS_6m.y
  
  # put block back toghether
  auto<-rbind(level4, level3, level2, level1)
  auto<-auto[order(auto$id, auto$side, -auto$difnum),]
  
  # keep only data needed for further steps
  below<-subset(auto, select=c(id, side, Md, Motor_2w, MS_2w, MS_6m_auto, MS_6m))
  
  # format data for further steps
  welow<-below
  welow$id <- factor(below$id)
  welow$side <- factor(below$side)
  welow$Md <- droplevels(below$Md)
  welow$MS_2w <- factor(below$MS_2w, levels=0:5, ordered=TRUE)
  welow$MS_6m <- factor(below$MS_6m, levels=0:5, ordered=TRUE)
  welow$MS_6m_auto <- factor(below$MS_6m_auto, levels=0:5, ordered=TRUE)
  
  # set motor scores model contrasts
  welow$Md<-relevel(welow$Md, ref="C5.-1")
  contrasts(welow$MS_2w)<-contr.treatment(6, base=1, contrasts=TRUE)
  contrasts(welow$MS_6m_auto)<-contr.treatment(6, base=1, contrasts=TRUE)
  
  # fit reference model for simulation
  welow$treat_num<-rnorm(dim(welow)[1],0,0.1)
  reference_polr_model<-polr(MS_6m ~ Md + MS_2w + MS_6m_auto + treat_num, data=welow, Hess=TRUE)

  return(reference_polr_model)}


### generate baseline data #################################################

generate_2w<-function(total_trialsize=50){
  
  # select only distance 0 (at Motor Level) for further manipulation
  man<-subset(mydata, diff=="0")
  
  # extract Motor Level left.right
  motl<-apply(  tapply(man$Motor_2w, list(man$id, man$side), paste), 1, function(x) paste0(x[1],sep=".", x[2]))
  
  # proportion of Motor Level combiantions
  pp<-prop.table(table(motl))
  # set T1.T1 combination to zero, as it does not contribute to the analysis!
  pp["T1.T1"]<-0
  
  # generate motor Level left.right 
  simpatt<-rmultinom(factor(names(pp)), size=total_trialsize, prob=pp)
  pat<-rep( factor(names(pp)),  simpatt)
  
  # sample EMSCI patients with given Motor Levels
  trial_baseline<-do.call( rbind.data.frame, 
                          lapply(pat, function(x) (mydata[mydata$id==sample(names(motl[motl==x]), size=1),]))  )
  
  # create simulation id (because same id may occurr multiple times)
  segs<- c(10,9,8,7,6,  9,8,7,6,5,  8,7,6,5,   7,6,5,4,3,   6,5,3,2 )
  trial_baseline$sim_id<-rep( 1:total_trialsize, rep(segs, simpatt ) )
  
  # decide trial arm (1:1 allocation rate)
  whotreated<-sample(unique(trial_baseline$sim_id), size=total_trialsize/2 )
  trial_baseline$treatment<-as.factor(ifelse(trial_baseline$sim_id %in% whotreated, "Y", "N"))
  
  # treatment as numeric for later steps
  trial_baseline$treat_num<-(1:nlevels(trial_baseline$treatment)-1)[trial_baseline$treatment]
  
  return(trial_baseline) }


### model based prediction of 6 months Motor Scores #####################################################

mob_preds<-function(tot_trialsize=10, treatment_effect=0){
  
  # draw Motor Level and sample patient for baseline data
  mybase<-generate_2w(total_trialsize=tot_trialsize)
  
  # Motor scores probbabilities at Motor Level 
  atlevel<-subset(mydata, difnum==0, keep= -diff)
  levprob<-prop.table(table(atlevel$MS_6m))
  
  # prediction of Motor Scores at 6 months for key muscle at Motor Level
  at<-subset(mybase, difnum==0)
  at$MS_6m_pred<- apply(at, 1, function(x) { t(0:5) %*% rmultinom(n=1, size=1, prob=levprob) } )
  
  # for subsequent steps you need to consider sides separately
  atr<-subset(at, difnum==0 & side == "right")
  atl<-subset(at, difnum==0 & side == "left")
  
  # laod reference model from "ref_mod.R"
  reference_polr_model<-ref_mod()
  
  # simulated treatment effect
  reference_polr_model$coefficients["treat_num"]<-treatment_effect 
  
  # model-based prediction of Motor scores for key muscles below the Motor Level
  draw_MS<-function(x) { t(0:5) %*% rmultinom(n=1, size=1, 
                        prob= predict(reference_polr_model, newdata=list("Md"=x[1], "MS_2w"=x[2], 
                        "MS_6m_auto"=x[4], "treat_num"=as.numeric(x[3])), type="prob") ) }
  
  # distance from Motor Level -1
  below1r<-subset(mybase, diff == -1 & side=="right")
  below1r$MS_6m_auto<-ordered(atr$MS_6m_pred[atr$sim_id %in% below1r$sim_id])
  below1r$MS_6m_pred<-apply(below1r[, c("Md", "MS_2w","treat_num", "MS_6m_auto")] ,1, draw_MS )
  
  below1l<-subset(mybase, diff == -1 & side=="left")
  below1l$MS_6m_auto<-ordered(atl$MS_6m_pred[atl$sim_id %in% below1l$sim_id])
  below1l$MS_6m_pred<-ordered( apply( below1l[, c("Md", "MS_2w","treat_num", "MS_6m_auto")] ,1, draw_MS ) , levels=0:5)
  
  # distance from Motor Level -2
  below2r<-subset(mybase, diff == -2 & side=="right")
  below2r$MS_6m_auto<-ordered(below1r$MS_6m_pred[below1r$sim_id %in% below2r$sim_id])
  below2r$MS_6m_pred<-ordered( apply( below2r[, c("Md", "MS_2w","treat_num", "MS_6m_auto")] ,1, draw_MS ) , levels=0:5)
  
  below2l<-subset(mybase, diff == -2 & side=="left")
  below2l$MS_6m_auto<-ordered(below1l$MS_6m_pred[below1l$sim_id %in% below2l$sim_id])
  below2l$MS_6m_pred<-ordered( apply( below2l[, c("Md", "MS_2w","treat_num", "MS_6m_auto")] ,1, draw_MS ) , levels=0:5)
  
  # distance from Motor Level -3
  below3r<-subset(mybase, diff == -3 & side=="right")
  below3r$MS_6m_auto<-ordered(below2r$MS_6m_pred[below2r$sim_id %in% below3r$sim_id])
  below3r$MS_6m_pred<-ordered( apply( below3r[, c("Md", "MS_2w","treat_num", "MS_6m_auto")] ,1, draw_MS ) , levels=0:5)
  
  below3l<-subset(mybase, diff == -3 & side=="left")
  below3l$MS_6m_auto<-ordered(below2l$MS_6m_pred[below2l$sim_id %in% below3l$sim_id])
  below3l$MS_6m_pred<-ordered( apply( below3l[, c("Md", "MS_2w","treat_num", "MS_6m_auto")] ,1, draw_MS ) , levels=0:5)
  
  # distance from Motor Level -4
  below4r<-subset(mybase, diff == -4 & side=="right")
  below4r$MS_6m_auto<-ordered(below3r$MS_6m_pred[below3r$sim_id %in% below4r$sim_id])
  below4r$MS_6m_pred<-ordered( apply( below4r[, c("Md", "MS_2w","treat_num", "MS_6m_auto")] ,1, draw_MS ) , levels=0:5)
  
  below4l<-subset(mybase, diff == -4 & side=="left")
  below4l$MS_6m_auto<-ordered(below3l$MS_6m_pred[below3l$sim_id %in% below4l$sim_id])
  below4l$MS_6m_pred<-ordered( apply( below4l[, c("Md", "MS_2w","treat_num", "MS_6m_auto")] ,1, draw_MS ) , levels=0:5)
  
  # data formatting  
  atr$MS_6m_auto<-NA
  atl$MS_6m_auto<-NA
  trial<-rbind(atr, atl, below1r,below1l,below2r,below2l,below3r,below3l,below4r,below4l)
  
  ttrial<-trial
  ttrial$id <- factor(trial$id)
  ttrial$side <- factor(trial$side)
  ttrial$Md <- factor(trial$Md)
  ttrial$treatment <- factor(trial$treatment)
  ttrial$diff <- droplevels(trial$diff)
  ttrial$sim_id <- factor(trial$sim_id)
  ttrial$Motor_2w <- factor(trial$Motor_2w)
  ttrial$MS_2w <- factor(trial$MS_2w, levels=0:5, ordered=TRUE)
  ttrial$MS_6m <- factor(trial$MS_6m, levels=0:5, ordered=TRUE)
  ttrial$MS_6m_pred <- factor(trial$MS_6m_pred, levels=0:5, ordered=TRUE)
  ttrial$MS_6m_auto <- factor(trial$MS_6m_auto, levels=0:5, ordered=TRUE)
  
  trial_data<-ttrial[order(ttrial$sim_id, ttrial$side, -as.numeric(ttrial$diff)),]
  
  return(trial_data)
  
}


### data generating proces  #####################################################

# create datasets
simdat <- lapply(1:nrow(param), function(i) { mob_preds(param[i, "n"], param[i, "treat"])})

# save datasets to simdat.rda
names(simdat) <- rownames(param)
save(param, simdat, file = "simdat.rda")


### tests battery  #####################################################

### t_test_6mos.R #############################################

t_test_6mos<-function( datain=my_complete_sim ){
  
  # select only relevant variables from simulated trial
  mydat<-subset(datain, select=c(id, sim_id, side, Motor_2w, MS_2w, MS_6m_pred, treatment))
  
  # convert Motor Scores to numeric variables
  mydat$MS_6m_pred<-(1 : nlevels(mydat$MS_6m_pred) -1)[mydat$MS_6m_pred]
  # do not compute difference total UEMS_6m - total UEMS_2w
  
  # sum all key segments
  MS_6mos<-cbind(sim_id=unique(mydat$sim_id), 
                 tot_MS_6mos=tapply(mydat$MS_6m_pred, mydat$sim_id, sum))
  
  # get other relevant information
  pat_info<-mydat[!duplicated(mydat$sim_id),c("sim_id", "treatment")]
  
  # merge 
  data_for_test<-merge(MS_6mos, pat_info, by="sim_id")
  
  # perform test and return p-value
  p<-tryCatch(t.test(tot_MS_6mos ~ treatment, data=data_for_test)$p.value, error=function(e) NA)
  
  return(as.vector(p))
}


### t_test_delta.R #############################################

t_test_delta<-function( datain=my_complete_sim ){
  
  # select only relevant variables from simulated trial
  mydat<-subset(datain, select=c(id, sim_id, side, Motor_2w, MS_2w, MS_6m_pred, treatment))
  
  # convert Motor Scores to numeric variables to compute difference
  mydat$MS_2w<-(1 : nlevels(mydat$MS_2w) -1)[mydat$MS_2w]
  mydat$MS_6m_pred<-(1 : nlevels(mydat$MS_6m_pred) -1)[mydat$MS_6m_pred]
  
  # compute difference total UEMS_6m - total UEMS_2w
  delta_MS<-cbind(sim_id=unique(mydat$sim_id), 
                  delta_MS=tapply(mydat$MS_6m_pred, mydat$sim_id, sum) - tapply(mydat$MS_2w, mydat$sim_id, sum))
  
  # get other relevant information
  pat_info<-mydat[!duplicated(mydat$sim_id),c("sim_id", "treatment")]
  
  # merge 
  data_for_test<-merge(delta_MS, pat_info, by="sim_id")
  
  # perform test and return p-value
  p<-tryCatch(t.test(delta_MS ~ treatment, data=data_for_test)$p.value, error=function(e) NA)
  
  return(as.vector(p))
}


### i_test_6mos.R #############################################

i_test_6mos<-function( datain=my_complete_sim ){
  
  # select only relevant variables from simulated trial
  mydat<-subset(datain, select=c(id, sim_id, side, Motor_2w, MS_2w, MS_6m_pred, treatment))
  
  # convert Motor Scores variables to compute operations
  mydat$MS_6m_pred<-(1 : nlevels(mydat$MS_6m_pred) -1)[mydat$MS_6m_pred]
  
  # compute total UEMS_6m
  MS<-cbind(sim_id=unique(mydat$sim_id),
            total_MS=tapply(mydat$MS_6m_pred, mydat$sim_id, sum))
  
  # get other relevant information
  pat_info<-mydat[!duplicated(mydat$sim_id),c("sim_id", "Motor_2w", "treatment")]
  # Motor_2w: if unsymmetric, it takes left, as it comes first into data
  
  # merge 
  data_for_test<-merge(MS, pat_info, by="sim_id")
  
  # excluded Motor Level for which there is not at least 1 patients for each treatment
  data_for_test_adj<-subset(data_for_test, 
                            Motor_2w %in% names( which( apply(table(Motor_2w, treatment), 1, min) > 1) )) 
  
  data_for_test_adj$Motor_2w<-droplevels(data_for_test_adj$Motor_2w)
  
  # perform independence test and return p-value
  itp<-tryCatch(pvalue(independence_test(  total_MS ~ treatment | Motor_2w , data=data_for_test_adj))[1], error=function(e) NA)
  
  return(itp)
}


### i_test_delta.R #############################################

i_test_delta<-function( datain=my_complete_sim ){
  
  # select only relevant variables from simulated trial
  mydat<-subset(datain, select=c(id, sim_id, side, Motor_2w, MS_2w, MS_6m_pred, treatment))
  
  # convert Motor Scores variables to compute operations
  mydat$MS_2w<-(1 : nlevels(mydat$MS_2w) -1)[mydat$MS_2w]
  mydat$MS_6m_pred<-(1 : nlevels(mydat$MS_6m_pred) -1)[mydat$MS_6m_pred]
  
  # compute difference total UEMS_6m - total UEMS_2w
  delta_MS<-cbind(sim_id=unique(mydat$sim_id),
                  delta_MS=tapply(mydat$MS_6m_pred, mydat$sim_id, sum) - tapply(mydat$MS_2w, mydat$sim_id, sum))
  
  # get other relevant information
  pat_info<-mydat[!duplicated(mydat$sim_id),c("sim_id", "Motor_2w", "treatment")]

  # merge 
  data_for_test<-merge(delta_MS, pat_info, by="sim_id")
  
  # excluded Motor Level for which there is not at least 1 patients for each treatment
  data_for_test_adj<-subset(data_for_test, Motor_2w %in% names( which( apply(table(Motor_2w, treatment), 1, min) > 1) )) 
  
  data_for_test_adj$Motor_2w<-droplevels(data_for_test_adj$Motor_2w)
  
  # perform independence test and return p-value
  itp<-tryCatch(pvalue(independence_test(  delta_MS ~ treatment | Motor_2w , data=data_for_test_adj))[1], error=function(e) NA)
  
  return(itp)
}


### a_test.R #############################################

a_test<-function( datain=my_complete_sim ){
  
  # select only relevant variables from simulated trial
  mydat<-subset(datain, select=c(id, sim_id, side, Motor_2w, MS_2w, MS_6m_pred, treatment))
  
  # convert Motor Scores to numeric variables to compute difference
  mydat$MS_2w<-(1 : nlevels(mydat$MS_2w) -1)[mydat$MS_2w]
  mydat$MS_6m_pred<-(1 : nlevels(mydat$MS_6m_pred) -1)[mydat$MS_6m_pred]
  
  # compute total UEMS_6m, and total UEMS_2w
  totals_MS<-cbind(sim_id=unique(mydat$sim_id), 
                   total_MS_6mos=tapply(mydat$MS_6m_pred, mydat$sim_id, sum),
                   total_MS_2weeks=tapply(mydat$MS_2w, mydat$sim_id, sum))
  
  # get other relevant information
  pat_info<-mydat[!duplicated(mydat$sim_id),c("sim_id", "treatment")]
  
  # merge 
  data_for_test<-merge(totals_MS, pat_info, by="sim_id")
  
  # perform test and return p-value
  p<-tryCatch( coef(summary(lm(total_MS_6mos ~ total_MS_2weeks + treatment, data=data_for_test)))[[3,4]], error=function(e) NA)

  return(as.vector(p))
}


### p_test.R #############################################

# define stripped POLR function for quicker simulation 
POLR <- function(m) {             
  
  x <- model.matrix(m)[,-1]
  x <- x[, colnames(x) %in% names(coef(m))]
  y <- unclass(model.response(model.frame(m)))
  
  s <- c(coef(m)[-length(coef(m))], 0,  m$zeta)
  names(s)[length(coef(m))]<-names(coef(m)[length(coef(m))])
  
  o <- rep(0, nrow(x))   
  me <- "logistic"
  
  function(which, perm) {
    x[,which] <- as.numeric(perm)
    MASS:::polr.fit(x, y, rep(1, nrow(x)), s, o, me, hessian = FALSE)$res$par[which]
  }
}

# define test

p_test<-function( datain=my_complete_sim, niter=1000 ){
  
  # consider segments below the motor level
  din<-subset(datain, difnum < 0)
  
  # set contrasts
  contrasts(din$MS_2w)<-contr.treatment(6, base=1, contrasts=TRUE)
  contrasts(din$MS_6m_auto)<-contr.treatment(6, base=1, contrasts=TRUE)
  
  # relevel Md
  din$Md<-relevel(din$Md, ref="C5.-1")
  
  # fit model
  model <- tryCatch(polr(MS_6m_pred ~ Md + MS_2w + MS_6m_auto + treatment, data=din, Hess=FALSE),
    error=function(e) "fitting_error")
  
  # error handling
  if( class(model) == "polr" )
    
  {
    obs_polr_beta <-coef(model)["treatmentY"]
    
    pfun <- POLR(model)
    
    xt <- (table(datain$sim_id, datain$treatment)[,"Y"] > 0)
    s <- names(xt)
    
    betas <- as.numeric(replicate(niter, tryCatch(
      pfun("treatmentY", din$sim_id %in% s[sample(xt)] ),
      error=function(e) NA )))
    
    # return p-value for observed treatment effect
    p<-mean(abs(betas) >= abs(obs_polr_beta), na.rm = TRUE )
    return(as.vector(p))}
  
  else 
  {return(as.vector(p<-NA))} 
}


### run tests on simulated data #############################################

# t-test on UEMS at 6 months
p_t_6mos <- unlist(lapply(simdat, t_test_6mos))

# t-test on UEMS delta change
p_t_delta <- unlist(lapply(simdat, t_test_delta))

# blocked independence test on UEMS 6 motnhs
p_i_6mos <- unlist(lapply(simdat, i_test_6mos))

# blocked independence test on delta change UEMS
p_i_delta <- unlist(lapply(simdat, i_test_delta))

# ANCOVA test
p_a <- unlist(lapply(simdat, a_test))

# transitional ordinal test
p_p <- unlist(lapply(simdat, p_test)) 


### Organize and store results #############################################

param$n <- ordered(param$n)
param$treat <- ordered(param$treat)

# import p-values from tests
param$t_test_6mos <- p_t_6mos
param$t_test_delta <- p_t_delta
param$i_test_6mos <- p_i_6mos
param$i_test_delta <- p_i_delta
param$a_test <- p_a
param$p_test <- p_p

results_all<-param

# calculate % of rejected tests
settings<-levels(unique(interaction( param$n, param$treat, sep= ":")))

t_6mos_rej <-with(param, tapply(p_t_6mos < .05, interaction(n, treat, sep = ":"), function(x) mean(x, na.rm = TRUE)))
t_delta_rej <-with(param, tapply(p_t_delta < .05, interaction(n, treat, sep = ":"), function(x) mean(x, na.rm = TRUE)))
i_6mos_rej <-with(param, tapply(p_i_6mos < .05, interaction(n, treat, sep = ":"), function(x) mean(x, na.rm = TRUE)))
i_delta_rej <-with(param, tapply(p_i_delta < .05, interaction(n, treat, sep = ":"), function(x) mean(x, na.rm = TRUE)))
a_rej <-with(param, tapply(p_a < .05, interaction(n, treat, sep = ":"), function(x) mean(x, na.rm = TRUE)))
p_rej <-with(param, tapply(p_p < .05, interaction(n, treat, sep = ":"), function(x) mean(x, na.rm = TRUE)))

# store power for each settings combination
results_summary<-as.data.frame(cbind(settings, t_6mos_rej, t_delta_rej, i_6mos_rej,  i_delta_rej, a_rej, p_rej ))

# save simulation results
save(results_all, results_summary, file="sim_results.rda")




