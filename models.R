# Construct time-only numerator and covariate-adjusted denominator models.
load_data <- function(prepared) {
  panel <- copy(prepared$panel)
  geom <- copy(prepared$geom)
  xx <- copy(prepared$xx)
  clones <- copy(prepared$clones)
  clones$artificial_censor <- tolower(as.character(clones$artificial_censor)) %in% c('true','1')
  people <- sort(unique(clones$stay_id)); N <- length(people)
  for(z in c('panel','geom','clones')) set(get(z), j='person', value=match(get(z)$stay_id,people))
  geom[,A:=as.integer(arm=='A25')]
  gr <- match(geom$row_id,panel$row_id)
  lastidx <- which(panel$last==1); lastidx <- lastidx[match(clones$clone,panel$clone[lastidx])]
  stopifnot(nrow(clones)==2*N,!anyNA(gr),!anyNA(lastidx),sum(geom$event)==sum(clones$analysis_event))
  features <- setdiff(names(xx),'covariate_id')
  cov <- as.matrix(xx[match(panel$covariate_id,xx$covariate_id),..features])
  stopifnot(all(is.finite(cov)))
  models <- list()
  for(phase in unique(panel$model)) {
    ii <- which(panel$model==phase); z <- panel[ii]
    tb <- if(phase=='B_wait') ns(log1p(z$clock_h),df=3) else matrix(z$clock_h,ncol=1)
    colnames(tb) <- paste0('clock_ns',seq_len(ncol(tb)))
    dose_event <- z$censor_next==1 & z$end_reason %in% c('NE_ge035_before_VP','NE_gt050_before_VP','VP_outside_dose_band')
    second <- if(phase=='B_wait') 'early_VP' else 'deadline'
    for(mechanism in c('dose',second)) for(kind in c('numerator','denominator')) {
      event <- if(mechanism=='dose') dose_event else z$censor_next==1 & z$end_reason==if(phase=='B_wait') 'VP_before_trigger' else 'no_VP_by_deadline'
      eligible <- z$fit_eligible==1 & (mechanism=='dose' | !dose_event)
      x <- cbind('(Intercept)'=1,tb,if(kind=='denominator') cov[ii,,drop=FALSE])
      center <- colMeans(x); scale <- apply(x,2,sd); center[1] <- 0
      scale[!is.finite(scale)|scale<1e-12] <- 1
      x <- sweep(sweep(x,2,center),2,scale,'/')
      models[[paste(phase,mechanism,kind,sep='_')]] <- list(ii=ii,x=x,use=eligible,y=as.integer(!event),person=z$person,kind=kind)
    }
  }
  stopifnot(length(models)==12)
  list2env(mget(c('N','people','panel','geom','xx','clones','cov','features','models','gr','lastidx'),envir=environment()),envir=.GlobalEnv)
  invisible(N)
}
build_km_geometry <- function() {
lapply(setNames(c('A25','B35'),c('A25','B35')),function(a){
 ii<-which(geom$arm==a);z<-geom[ii];tt<-sort(unique(c(z$tstart,z$tstop)))
 dd<-sort(unique(z$tstop[z$event==1]));nd<-length(dd)
 change<-sparseMatrix(i=c(match(z$tstart,tt),match(z$tstop,tt)),j=rep(seq_len(nrow(z)),2),
   x=c(rep(1,nrow(z)),rep(-1,nrow(z))),dims=c(length(tt),nrow(z)))
 event<-sparseMatrix(i=match(z$tstop[z$event==1],dd),j=which(z$event==1),x=1,dims=c(nd,nrow(z)))
 list(ii=ii,change=change,event=event,death_times=dd,before=findInterval(dd,tt,left.open=TRUE))
})

}
