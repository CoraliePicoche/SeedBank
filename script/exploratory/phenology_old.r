
start_duration_bloom=function(time_series_inter,date_bis,hydro_interp,threshold){ #This function is supposed to detect the date of the beginning and duration of bloom, defined as the population going above a certain threshold/quantile, with a gradient positive (beginning) or negative.
#	dates_bis=seq(dates[1],dates[length(dates)],1) #Regular time grid
#	time_series_interp=na.approx(time_series,x=dates,xout=dates_bis,na.rm=FALSE)
	gradient=diff(time_series_interp)
	ab_quantile=quantile(time_series_interp,threshold,na.rm=T)
	dates_above=which(time_series_interp>ab_quantile)
	beg=c()
	end_tmp=c()
#	if(dates_above[length(dates_above)]<length(dates_bis)){
#		end_tmp=c(dates_bis[dates_above[length(dates_above)]+1])
#	}else{
#		end_temp=dates_bis[length(dates_bis)]+1
#	}
	for(i in 1:length(dates_above)){
		if(dates_above[i]==1){
			beg=c(beg,dates_bis[dates_above[i]])
		}else{
			if(!is.na(gradient[dates_above[i]-1])){
				if((dates_bis[dates_above[i]]+1)<dates_bis[length(dates_bis)]){
				if((gradient[dates_above[i]-1]>0)&(time_series_interp[dates_above[i]-1]<ab_quantile)&((time_series_interp[dates_above[i]+1]>ab_quantile))){
					beg=c(beg,dates_bis[dates_above[i]])
				}
				else if ((gradient[dates_above[i]-1]<0)&(time_series_interp[dates_above[i]-1]>ab_quantile)&(time_series_interp[dates_above[i]+1]<ab_quantile)){
#					if(end_tmp[length(end_tmp)]+1!=dates_bis[dates_above[i]]){
						end_tmp=c(end_tmp,dates_bis[dates_above[i]])
#					}
	
				}
			}else{
                                if((gradient[dates_above[i]-1]>0)&(time_series_interp[dates_above[i]-1]<ab_quantile)){
                                        beg=c(beg,dates_bis[dates_above[i]])
                                }
                                else if ((gradient[dates_above[i]-1]<0)&(time_series_interp[dates_above[i]-1]>ab_quantile)){
                                                end_tmp=c(end_tmp,dates_bis[dates_above[i]])

                                }

			}
		}
		}
	}
	
	dd_tmp=dates_bis[!is.na(time_series_inter)]
	dd1=dd_tmp[1] #First date for which we have no NA
	dd2=dd_tmp[length(dd_tmp)] #Final date for which we have no NA

	#Now compute duration
	only_beg=beg[1]
	only_end=end_tmp[length(end_tmp)]
	temp_beg=c()
        duration=c()
        if(length(beg)==2){
                if(beg[2]-beg[1]>2){
                        if((as.Date(beg[1])!=dd1)&&(as.Date(beg[2])!=dd2)){
                                duration=c(duration,beg[i]-only_beg[length(only_beg)])
                        }
                        only_beg=c(only_beg,beg[2])
                }
        }else{
                for(i in 2:(length(beg)-1)){
                if((beg[i+1]-beg[i])>2){
                        if((as.Date(beg[i])!=dd1)&&(as.Date(beg[i+1])!=dd2)){
                                duration=c(duration,beg[i]-only_beg[length(only_beg)])
                        }
                        only_beg=c(only_beg,beg[i+1])
                #       only_end=
                        temp_beg=c(temp_beg,hydro_interp[dates_bis==beg[i]])
                }
        }
        }
        temp_beg=c(temp_beg,hydro_interp[dates_bis==beg[i]])

        return(list(only_beg,duration,end_tmp,temp_beg))
}

