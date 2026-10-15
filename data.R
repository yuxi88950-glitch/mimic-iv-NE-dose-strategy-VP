# Independently simulated patients; no clinical records or fitted study estimates.
generate_data <- function(n=800L, seed=20261016L) {
  set.seed(seed)
  severity <- rnorm(n)
  people <- data.table(stay_id=seq_len(n),
    age=pmin(95,pmax(18,round(rnorm(n,65,14)))), male=rbinom(n,1,.55),
    baseline_sofa=pmin(22,pmax(2,round(10+2*severity+rnorm(n,0,2)))),
    ne0=runif(n,.25,.349), episode_end=runif(n,18,48))
  death <- rexp(n,exp(.4*severity)/850)
  people[, `:=`(followup_h=pmin(death,672),event28=as.integer(death<=672),
                first_vp=Inf,documented_stop=Inf)]
  observations <- vector('list',n)
  for(i in seq_len(n)) {
    p <- people[i]; time <- seq(0,floor(p$episode_end*2)/2,by=.5)
    ne <- numeric(length(time)); ne[1] <- p$ne0
    drift <- sample(c(-.012,.006,.016),1,prob=c(.35,.35,.30))+.002*severity[i]
    for(k in 2:length(time)) {
      if(ne[k-1]==0) ne[k] <- 0 else {
        ne[k] <- max(0,ne[k-1]+drift+rnorm(1,0,.017)+ifelse(runif(1)<.02,.24,0))
        if(ne[k]<.03) ne[k] <- 0
      }
    }
    map <- 70-4*severity[i]-12*(ne-.3)+rnorm(length(time),0,5)
    lactate <- exp(1+.3*severity[i]+rnorm(length(time),0,.15))
    creatinine <- exp(.1+.35*severity[i]+rnorm(length(time),0,.10))
    lactate[time%%6!=0 | runif(length(time))<.2] <- NA_real_
    creatinine[time%%12!=0 | runif(length(time))<.1] <- NA_real_
    sofa <- pmin(24,pmax(0,p$baseline_sofa+round(5*(ne-p$ne0))))
    ventilated <- rbinom(1,1,plogis(.3+.5*severity[i]))
    other <- as.integer(ne>.4 & severity[i]>0)
    candidates <- which(time>0 & ne>0 & time<p$followup_h &
      runif(length(time))<plogis(-2.8+7*(ne-.3)+.3*severity[i]+.4*(map<65)))
    if(length(candidates)) people$first_vp[i] <- time[candidates[1]]
    off <- which(ne==0)
    if(length(off)) people$documented_stop[i] <- time[off[1]]
    observations[[i]] <- data.table(stay_id=i,time=time,ne=ne,map=map,
      lactate=lactate,creatinine=creatinine,sofa=sofa,
      ventilated=ventilated,other_vasoactive=other)
  }
  list(people=people,observations=rbindlist(observations))
}

# Reconstruct the two strategies using only the observed synthetic history.
prepare_data <- function(data) {
  people <- data$people
  histories <- split(data$observations,data$observations$stay_id)
  panel_list <- geometry_list <- clone_list <- matrix_list <- list()
  nr <- nc <- nx <- ng <- 0L
  band <- function(x) ifelse(x<.25,'below',ifelse(x<.35,'low',ifelse(x<=.5,'high','above')))
  first_time <- function(z,condition) if(any(condition)) min(z$time[condition]) else Inf
  covariates <- function(p,z,t) {
    j <- max(which(z$time<=t)); current <- z[j]
    value <- function(name,hours) {
      use <- which(z$time<=t & z$time>=t-hours & !is.na(z[[name]]))
      if(length(use)) z[[name]][tail(use,1)] else NA_real_
    }
    previous <- z$ne[max(which(z$time<=max(0,t-1)))]
    x <- data.table(age=p$age,male=p$male,baseline_sofa=p$baseline_sofa,
      map=value('map',2),lactate=value('lactate',6),creatinine=value('creatinine',24),
      ne=current$ne,sofa=value('sofa',24),ventilated=current$ventilated,
      other_vasoactive=current$other_vasoactive,
      ne_from_t0=current$ne-p$ne0,ne_change_1h=current$ne-previous)
    for(v in c('map','lactate','creatinine','sofa')) {
      x[[paste0(v,'_missing')]] <- as.integer(is.na(x[[v]]))
      x[[v]][is.na(x[[v]])] <- 0
    }
    x
  }
  for(i in seq_len(nrow(people))) {
    p <- people[i]; z <- histories[[as.character(p$stay_id)]]
    stopifnot(p$ne0>=.25,p$ne0<.35)
    for(arm in c('A25','B35')) {
      nc <- nc+1L; clone <- nc
      trigger <- if(arm=='A25') 0 else first_time(z,z$ne>=.35 & z$ne<=.5)
      upper <- first_time(z,if(arm=='A25') z$ne>=.35 else z$ne>.5)
      deadline <- trigger+3
      ev <- data.table(time=c(p$followup_h,p$documented_stop,p$episode_end,upper,deadline),
        priority=c(0,30,40,20,50),
        reason=c('death_before_deviation','documented_stop','episode_end',
                 if(arm=='A25') 'NE_ge035_before_VP' else 'NE_gt050_before_VP','no_VP_by_deadline'),
        kind=c('natural','natural','natural','artificial','artificial'))
      # A surviving patient has administrative follow-up, not a death at day 28.
      if(!p$event28) ev <- ev[reason!='death_before_deviation']
      if(is.finite(p$first_vp)) {
        vp <- p$first_vp
        # The NE rate effective at this timestamp is applied before VP classification.
        vp_band <- band(z$ne[max(which(z$time<=vp))])
        if(arm=='B35' && vp_band=='above') reason <- 'VP_outside_dose_band'
        else if(arm=='B35' && vp<trigger) reason <- 'VP_before_trigger'
        else if(vp_band!=if(arm=='A25') 'low' else 'high') reason <- 'VP_outside_dose_band'
        else if(vp>=trigger && vp<=deadline) reason <- 'VP_in_window_and_band'
        else reason <- 'VP_after_deadline'
        ev <- rbind(ev,data.table(time=vp,priority=10,reason=reason,
          kind=if(reason=='VP_in_window_and_band') 'completed' else 'artificial'))
      }
      setorder(ev,time,priority); terminal <- ev[1]
      end <- terminal$time; artificial <- terminal$kind=='artificial'
      if(trigger>end) trigger <- Inf
      analysis_stop <- if(artificial) min(end,p$followup_h) else p$followup_h
      analysis_event <- as.integer(p$event28==1 && p$followup_h<=analysis_stop)
      clone_list[[nc]] <- data.table(clone=clone,stay_id=p$stay_id,arm=arm,
        trigger_h=trigger,history_end_h=end,terminal_reason=terminal$reason,
        terminal_type=terminal$kind,artificial_censor=artificial,
        analysis_stop_h=analysis_stop,analysis_event=analysis_event)
      # Waiting uses T0; each grace window uses its own trigger-reset hourly grid.
      wait_end <- min(end,trigger)
      cuts <- c(0,end,seq(0,wait_end,by=1))
      if(is.finite(trigger)) cuts <- c(cuts,(trigger+0:3)[trigger+0:3<=end])
      cuts <- sort(unique(cuts))
      instant <- terminal$kind=='completed' && trigger==end
      first_row <- nr+1L
      for(k in seq_len(length(cuts)-1L)) {
        left <- cuts[k]; right <- cuts[k+1L]
        waiting <- arm=='B35' && left<trigger
        ending <- right==end
        natural <- ending && terminal$kind=='natural'
        censored <- ending && artificial
        nr <- nr+1L; nx <- nx+1L
        matrix_list[[nx]] <- cbind(data.table(covariate_id=nx),covariates(p,z,left))
        panel_list[[nr]] <- data.table(row_id=nr,clone=clone,stay_id=p$stay_id,arm=arm,
          model=if(waiting) 'B_wait' else paste0(arm,'_grace'),start_h=left,end_h=right,
          clock_h=if(waiting) left else left-trigger,last=as.integer(ending && !instant),
          censor_next=as.integer(censored),fit_eligible=as.integer(!natural),
          apply_update=as.integer(!natural && !censored),
          vp_start_in_interval=as.integer(ending && terminal$kind=='completed' && !instant),
          instant_start=0L,end_reason=if(ending && !instant) terminal$reason else
            if(waiting && right==trigger) 'trigger' else 'interval',covariate_id=nx)
      }
      if(instant || end==0) {
        # Include immediate initiation without adding survival time or post-VP covariates.
        nr <- nr+1L
        if(arm=='B35' && end>0) row <- copy(panel_list[[nr-1L]]) else {
          nx <- nx+1L
          matrix_list[[nx]] <- cbind(data.table(covariate_id=nx),covariates(p,z,0))
          row <- data.table(row_id=nr,clone=clone,stay_id=p$stay_id,arm=arm,
            model=if(arm=='B35') 'B_wait' else 'A25_grace',start_h=0,end_h=0,clock_h=0,last=1L,
            censor_next=as.integer(artificial),fit_eligible=as.integer(terminal$kind!='natural'),
            apply_update=as.integer(terminal$kind=='completed'),
            vp_start_in_interval=as.integer(instant),instant_start=as.integer(instant),
            end_reason=terminal$reason,covariate_id=nx)
        }
        if(instant) row[, `:=`(row_id=nr,model=paste0(arm,'_grace'),start_h=end,end_h=end,clock_h=0,
                    last=1L,vp_start_in_interval=1L,instant_start=1L,
                    end_reason=terminal$reason)]
        panel_list[[nr]] <- row
      }
      for(j in first_row:nr) {
        r <- panel_list[[j]]; stop <- min(r$end_h,analysis_stop)
        if(stop>r$start_h) {
          ng <- ng+1L
          geometry_list[[ng]] <- data.table(clone=clone,stay_id=p$stay_id,arm=arm,
            row_id=r$row_id,tstart=r$start_h,tstop=stop,
            event=as.integer(analysis_event==1 && stop==analysis_stop),carry=0L)
        }
        if(r$last && analysis_stop>r$end_h) {
          ng <- ng+1L
          geometry_list[[ng]] <- data.table(clone=clone,stay_id=p$stay_id,arm=arm,
            row_id=r$row_id,tstart=r$end_h,tstop=analysis_stop,event=analysis_event,carry=1L)
        }
      }
    }
  }
  list(panel=rbindlist(panel_list),geom=rbindlist(geometry_list),
       xx=rbindlist(matrix_list),clones=rbindlist(clone_list))
}
