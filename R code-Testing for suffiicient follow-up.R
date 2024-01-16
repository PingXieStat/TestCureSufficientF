
####### The distribution of the failure time for susceptibles is a Weibull distribution

######## generate.data functin is used to generate the survival data
generate.data<-function(s.size,p,shape.p,scale.p,cpara,null.in){
  # s.size : the sample size
  # p : the noncure rate
  # shape.p : the shape parameter of Weibull distribution for failure distribution
  # scale.p : the scale parameter of Weibull distribution for failure distribution
  # cpara : the parameter of the censoring distribution
  # null.in : hypo.in indicates whether the null hypothesis is true or not
  
  if(null.in=="T"){
    condi<-TRUE
    while(condi){
      cure.indi<-I(runif(s.size,0,1)<=p)*1 # 1 indicates the subject is non-cured.
      T<-rep(0,s.size)
      T[cure.indi==0]<-Inf
      T[cure.indi==1]<-rweibull(length(cure.indi[cure.indi==1]),shape=shape.p,scale=scale.p)
      C<-rweibull(s.size,shape=1,scale=cpara)   # Med(C)/Med(T|T<Inf)(Median) = cpara/scale.p
      Y<-pmin(T,C)
      delta<-I(T<=C)*1
      data<-data.frame(Y,delta)
      ifelse(delta[which.max(Y)]==1,condi<-TRUE,condi<-FALSE)  
    }
    return(data)
  }else{
    condi<-TRUE
    while(condi){
      cure.indi<-I(runif(s.size,0,1)<=p)*1 # 1 indicates the subject is non-cured.
      T<-rep(0,s.size)
      T[cure.indi==0]<-Inf
      T[cure.indi==1]<-rweibull(length(cure.indi[cure.indi==1]),shape=shape.p,scale=scale.p)
      C<-runif(s.size,0,cpara)
      Y<-pmin(T,C)
      delta<-I(T<=C)*1
      data<-data.frame(Y,delta)
      ifelse(delta[which.max(Y)]==1,condi<-TRUE,condi<-FALSE)
    }
    return(data)
  }
  
}

### F0hat function is the KM estimator for the failure time distribution
F0hat<-function(data,t){
  # data : dataset
  # Y: the observed time
  # delta: the censoring indicator
  
  if(t < min(data$Y)){
    return(0)
  }else{
    
    data <- data[order(data$Y),]
    s.size <- length(data$Y)
    
    v <- (1 - data$delta/(s.size-1:s.size+1))
    out <- 1-prod(v[data$Y<=t])
    
    return(out)
  }
}

####### test.mine function is used to calculate the proposed test statistic Tn=p_G-p_n
test.mine<-function(data,epsilon){
  #data: dataset
  #ep: the optimal epsilon
  
  ep<-epsilon
  
  t.max<-max(data$Y); tf.max<-max(data$Y[which(data$delta==1)]) 
  
  F0<-F0hat(data,t.max)
  pGhat<-function(t){
    out<-F0hat(data,t.max-t)+(F0hat(data,t.max-t/2)-F0hat(data,t.max-t))^2/
      (2*F0hat(data,t.max-t/2)-F0hat(data,t.max-t)-F0hat(data,t.max))
    return(out)
  }
  
  pGhat.value<-pGhat(ep)
  
  if((ep<=t.max-tf.max)|(ep>t.max)){
    out<-NULL
  }else{
    if(is.na(pGhat.value)){
      out<-0
    }else{
      if(pGhat.value>1){
        out<-1-F0
      }else{
        out<-ifelse((pGhat.value<F0),0,pGhat.value-F0)
      }
    }
  }
  
  
  return(out)
}

###boots funcition is used to calculate the bootstrap test statistics. 
boots<-function(data,epsilon){
  #data: dataset
  #ep: the optimal epsilon
  data1<-data
  
  ep<-epsilon
  
  s.size<-dim(data)[1]
  
  
  training<-NULL;condi<-TRUE;ite<-1;
  while(condi){
    condit<-TRUE
    while(condit){
      
      #generating the bootstrap sample
      iid<-sample(1:s.size,s.size,replace = TRUE)
      data.B<-data[iid,]
      
      condit<-ifelse(data.B$delta[which.max(data.B$Y)]==1,TRUE,FALSE)
    }
    
    ## test.sta.boots: calculating the bootstrap test statistics
    test.sta.boots<-test.mine(data.B,epsilon=ep)
    
    training<-c(training,test.sta.boots)
    ite<-ite+1
    condi<-ifelse(ite==1000|length(training)==500,FALSE,TRUE) 
  }
  
  training1<-training-test.mine(data=data1,epsilon=ep)
  
  if(length(training1)<5){ return(NULL) }else{
    return(quantile(training1,0.95)) # the output is 95% quantile
  }
  
}


########## main.fun is used to calculate the type I error or power
main.fun<-function(i){
  
  s<-i+1
  out<-NULL
  
  ######################### 1 
  ############################################################################
  
  data.p1<-generate.data(s.size=s.siz[1],p=pp,shape.p=shape.pp,scale.p=scale.pp,cpara=cparaa,null.in=null.inn)
  s.size<-dim(data.p1)[1]
  
  #################### Our method
  
  t.max<-max(data.p1$Y);tf.max<-max(data.p1$Y[which(data.p1$delta==1)])
  
  epsi.opti<-ifelse(2*(t.max-tf.max)>=t.max,t.max,2*(t.max-tf.max)+rate1*(t.max-2*(t.max-tf.max)))
  
  t.mine<-test.mine(data=data.p1,epsilon=epsi.opti) # t.mine is the test statistic
  Crv<-boots(data=data.p1,epsilon=epsi.opti) # Crv is the critical value
  power<-I(t.mine>Crv)*1 # power is the type I error or power
  
  out<-c(out,power)
  
  ######################### 2 
  ############################################################################
  
  data.p1<-generate.data(s.size=s.siz[2],p=pp,shape.p=shape.pp,scale.p=scale.pp,cpara=cparaa,null.in=null.inn)
  s.size<-dim(data.p1)[1]
  
  #################### Our method
  
  t.max<-max(data.p1$Y);tf.max<-max(data.p1$Y[which(data.p1$delta==1)])
  
  epsi.opti<-ifelse(2*(t.max-tf.max)>=t.max,t.max,2*(t.max-tf.max)+rate1*(t.max-2*(t.max-tf.max)))
  
  t.mine<-test.mine(data=data.p1,epsilon=epsi.opti) # t.mine is the test statistic
  Crv<-boots(data=data.p1,epsilon=epsi.opti) # Crv is the critical value
  power<-I(t.mine>Crv)*1 # power is the type I error or power
  
  out<-c(out,power)
  
  ######################### 3
  ############################################################################
  
  data.p1<-generate.data(s.size=s.siz[3],p=pp,shape.p=shape.pp,scale.p=scale.pp,cpara=cparaa,null.in=null.inn)
  s.size<-dim(data.p1)[1]
  
  #################### Our method
  
  t.max<-max(data.p1$Y);tf.max<-max(data.p1$Y[which(data.p1$delta==1)])
  
  epsi.opti<-ifelse(2*(t.max-tf.max)>=t.max,t.max,2*(t.max-tf.max)+rate1*(t.max-2*(t.max-tf.max)))
  
  t.mine<-test.mine(data=data.p1,epsilon=epsi.opti) # t.mine is the test statistic
  Crv<-boots(data=data.p1,epsilon=epsi.opti) # Crv is the critical value
  power<-I(t.mine>Crv)*1 # power is the type I error or power
  
  out<-c(out,power)
  
  
  ######################### 4
  ############################################################################
  
  data.p1<-generate.data(s.size=s.siz[4],p=pp,shape.p=shape.pp,scale.p=scale.pp,cpara=cparaa,null.in=null.inn)
  s.size<-dim(data.p1)[1]
  
  #################### Our method
  
  t.max<-max(data.p1$Y);tf.max<-max(data.p1$Y[which(data.p1$delta==1)])
  
  epsi.opti<-ifelse(2*(t.max-tf.max)>=t.max,t.max,2*(t.max-tf.max)+rate1*(t.max-2*(t.max-tf.max)))
  
  t.mine<-test.mine(data=data.p1,epsilon=epsi.opti) # t.mine is the test statistic
  Crv<-boots(data=data.p1,epsilon=epsi.opti) # Crv is the critical value
  power<-I(t.mine>Crv)*1 # power is the type I error or power
  
  out<-c(out,power)
  
  
  if(length(out)==4){
    return(out)
  }else{
    return(NULL)
  }
}




################### An example

null.inn<-"T"
shape.pp<-1.5;scale.pp<-1.5;cparaa<-3
pp<-0.9; s.siz<-c(400,800,1200,1800)
rate1<-7/8

main.fun(1)
