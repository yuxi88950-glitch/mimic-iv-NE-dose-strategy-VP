# Stabilized IPC weights. Clinical coefficients alone receive ridge penalties.

lambda_grid <- 10^seq(3,-5,length.out=33)

fit_numerator<-function(m,frequency){
 f<-frequency[m$person];use<-m$use & f>0
 x<-m$x;q<-qr(x[use,,drop=FALSE],tol=1e-10,LAPACK=FALSE)
 keep<-sort(q$pivot[seq_len(q$rank)]);x<-x[,keep,drop=FALSE];warnings<-character()
 solve<-function(engine,maxit)withCallingHandlers(
  if(engine=='fastglm')fastglm(x[use,,drop=FALSE],m$y[use],family=binomial(),weights=as.numeric(f[use]),
    method=0L,tol=1e-8,maxit=maxit,firth=FALSE)
  else glm.fit(x[use,,drop=FALSE],m$y[use],family=binomial(),weights=as.numeric(f[use]),
    control=glm.control(maxit=maxit,epsilon=1e-8)),
  warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 good<-function(g)!is.null(g)&&isTRUE(g$converged)&&all(is.finite(g$coefficients))&&g$rank==ncol(x)
 engine<-'fastglm';fit<-tryCatch(solve(engine,150L),error=function(e)NULL)
 if(!good(fit)){engine<-'glm.fit';fit<-tryCatch(solve(engine,500L),error=function(e)NULL)}
 if(!good(fit))stop('Unpenalized logistic fit failed; no automatic penalty substitution')
 p<-plogis(as.numeric(x%*%fit$coefficients));stopifnot(all(is.finite(p)))
 list(p=pmin(pmax(p,1e-6),1-1e-6),raw=p,
  monitor=data.table(rows=sum(use),weighted_rows=sum(f[use]),censors=sum(f[use]*(1-m$y[use])),
    rank=ncol(x),original_columns=ncol(m$x),converged=fit$converged,engine=engine,
    clip_low=sum(p[use]<1e-6),clip_high=sum(p[use]>1-1e-6),
    max_abs_coefficient=max(abs(fit$coefficients)),warnings=paste(unique(warnings),collapse=' | ')),
  coefficients=data.table(term=colnames(x),coefficient=as.numeric(fit$coefficients)),
  removed=setdiff(colnames(m$x),colnames(x)))
}

fit_one<-function(m,frequency,folds){
 if(m$kind=='numerator'){
  a<-fit_numerator(m,frequency)
  a$monitor[,`:=`(lambda=NA_real_,lambda_at_boundary=FALSE,cv_deviance=NA_real_,jerr=0L)]
  a$cv<-NULL;return(a)
 }
 f<-frequency[m$person];use<-m$use & f>0
 x<-m$x[,colnames(m$x)!='(Intercept)',drop=FALSE]
 # Ridge retains linearly dependent clinical columns; only constant columns
 # lack information. glmnet internally omits those constants.
 penalty<-as.numeric(!grepl('^clock_ns',colnames(x)));warnings<-character()
 solve_cv<-function(family_spec,modified=FALSE)withCallingHandlers(cv.glmnet(x[use,,drop=FALSE],m$y[use],family=family_spec,
   weights=as.numeric(f[use]),alpha=0,lambda=lambda_grid,penalty.factor=penalty,
   foldid=folds[m$person[use]],standardize=TRUE,type.measure='deviance',grouped=TRUE,
   thresh=if(modified)1e-7 else 1e-8,maxit=200000,type.logistic=if(modified)'modified.Newton'else'Newton',parallel=FALSE),
   warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 fit<-tryCatch(solve_cv('binomial'),error=function(e){
   warnings<<-c(warnings,conditionMessage(e));NULL
 })
 incomplete<-function(w)any(grepl('Convergence for|error code|empty model|convergence issue',w,ignore.case=TRUE))
 if(is.null(fit) || incomplete(warnings)){
   full_failure<-is.null(fit)
   initial_warnings<-unique(warnings);warnings<-character()
   fit<-tryCatch(if(full_failure)solve_cv(binomial())else solve_cv('binomial',modified=TRUE),error=function(e){warnings<<-c(warnings,conditionMessage(e));NULL})
   if(is.null(fit) || incomplete(warnings)){
     retry_warnings<-warnings;warnings<-character()
     fit<-solve_cv(binomial())
     if(incomplete(warnings))stop('CV path did not converge after numerical retries')
     initial_warnings<-c(initial_warnings,retry_warnings)
   }
   warnings<-c(warnings,paste('RESOLVED CV path using same-objective numerical retry (modified Newton tolerance 1e-7 or generic IRLS):',paste(initial_warnings,collapse=' | ')))
 }
 stopifnot(fit$glmnet.fit$jerr==0,all(is.finite(fit$cvm)))
 p<-as.numeric(predict(fit,newx=x,s='lambda.min',type='response'))
 beta<-as.matrix(coef(fit,s='lambda.min'))[,1]
 stopifnot(all(is.finite(p)),all(is.finite(beta)))
 constant<-colnames(x)[apply(x[use,,drop=FALSE],2,function(v)max(v)-min(v)<1e-12)]
 list(p=pmin(pmax(p,1e-6),1-1e-6),raw=p,
   monitor=data.table(rows=sum(use),weighted_rows=sum(f[use]),censors=sum(f[use]*(1-m$y[use])),
     rank=NA_integer_,original_columns=ncol(m$x),converged=TRUE,engine='cv.glmnet ridge',
     clip_low=sum(p[use]<1e-6),clip_high=sum(p[use]>1-1e-6),
     max_abs_coefficient=max(abs(beta)),warnings=paste(unique(warnings),collapse=' | '),
     lambda=fit$lambda.min,lambda_at_boundary=fit$lambda.min%in%range(lambda_grid),
     cv_deviance=min(fit$cvm),jerr=fit$glmnet.fit$jerr),
   coefficients=data.table(term=names(beta),coefficient=as.numeric(beta)),removed=constant,
   cv=data.table(lambda=fit$lambda,cv_deviance=fit$cvm,cv_se=fit$cvsd))
}

make_folds<-function(frequency,seed){
 set.seed(seed);fold<-rep(NA_integer_,N)
 den<-models[grepl('denominator$',names(models))]
 sets<-lapply(den,function(m)unique(m$person[m$use & m$y==0]))
 sets<-c(sets,list(seq_len(N)))
 for(ids in sets){
  ids<-ids[frequency[ids]>0];todo<-ids[is.na(fold[ids])]
  if(length(todo))todo<-todo[sample.int(length(todo))]
  for(i in todo){
   counts<-vapply(1:5,function(k)sum(frequency[ids[!is.na(fold[ids]) & fold[ids]==k]]),numeric(1))
   best<-which(counts==min(counts));total<-vapply(best,function(k)sum(frequency[which(fold==k)]),numeric(1))
   fold[i]<-best[which.min(total)]
  }
 }
 stopifnot(!anyNA(fold[frequency>0]));fold
}


fit_all<-function(frequency,keep=FALSE,seed=2026091701L){
 folds<-make_folds(frequency,seed);cv<-list()
 pn<-pd<-rep(1,nrow(panel));mon<-coefs<-removed<-list()
 for(nm in names(models)){
  m<-models[[nm]];a<-tryCatch(fit_one(m,frequency,folds),
    error=function(e)stop(nm,': ',conditionMessage(e),call.=FALSE))
  if(endsWith(nm,'numerator'))pn[m$ii]<-pn[m$ii]*a$p else pd[m$ii]<-pd[m$ii]*a$p
  mon[[nm]]<-a$monitor[,model:=nm]
  if(keep){if(!is.null(a$cv))cv[[nm]]<-a$cv[,model:=nm];coefs[[nm]]<-a$coefficients[,model:=nm]
    removed[[nm]]<-data.table(model=rep(nm,length(a$removed)),term=a$removed)}
 }
 list(pn=pn,pd=pd,monitor=rbindlist(mon),coefficients=rbindlist(coefs),removed=rbindlist(removed),cv=rbindlist(cv),folds=folds)
}

make_weights<-function(fit,frequency=rep(1L,N)){
 inc<-ifelse(panel$apply_update==1 & frequency[panel$person]>0,log(fit$pn)-log(fit$pd),0)
 endlog<-ave(inc,panel$clone,FUN=cumsum);startlog<-endlog-inc
 stopifnot(all(is.finite(endlog)),all(abs(endlog)<700))
 sw<-exp(startlog);ew<-exp(endlog)
 final<-ifelse(clones$analysis_stop_h>panel$end_h[lastidx],ew[lastidx],sw[lastidx])
 retained<-!clones$artificial_censor
 sampled<-rep(final[retained],frequency[clones$person[retained]])
 limits<-unname(quantile(sampled,c(.01,.99),type=7));lower<-limits[1];upper<-limits[2]
 actual<-ifelse(geom$carry==1,ew[gr],sw[gr])
 list(weight=pmin(pmax(actual,lower),upper),raw=actual,lower=lower,upper=upper,final=final,start=sw,end=ew,increment=inc)
}

absolute_effects<-function(weights,frequency=rep(1,N),keep=FALSE){
 values<-numeric();curves<-list()
 for(a in names(km_geometry)){
  z<-km_geometry[[a]];mass<-weights[z$ii]*frequency[geom$person[z$ii]]
  cum<-cumsum(as.numeric(z$change%*%mass));risk<-c(0,cum)[z$before+1L]
  deaths<-as.numeric(z$event%*%mass);use<-deaths>0
  tt<-z$death_times[use];risk<-risk[use];deaths<-deaths[use]
  stopifnot(all(risk>0),all(deaths<=risk+1e-6))
  s<-cumprod(pmax(0,1-deaths/risk))
  values[paste0('risk_',a)]<-if(length(s))1-tail(s,1)else 0
  values[paste0('rmst_',a)]<-sum(diff(c(0,tt,672))*c(1,s))/24
  if(keep)curves[[a]]<-data.table(arm=a,time_h=c(0,tt,672),survival=c(1,s,if(length(s))tail(s,1)else 1))
 }
 values['RD']<-values['risk_A25']-values['risk_B35']
 values['RR']<-values['risk_A25']/values['risk_B35']
 values['RMST_difference']<-values['rmst_A25']-values['rmst_B35']
 list(values=values,curves=if(keep)rbindlist(curves)else NULL)
}

