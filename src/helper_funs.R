naif_nmin <- function(x, n_min) {
  if(length(na.omit(x)) > 0 && sum(!is.na(x)) >= n_min) {
    na.omit(x)
  } else NA
}

c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  # "black", remove black so can be used as reference
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

closest_point <- function(points, target) {
  if(is.vector(target)) {
    ind <- which.min(sp::spDistsN1(pts = as.matrix(points), pt = target))
    return(c(points[ind, 1], points[ind, 2]))
  }
  else if(is.data.frame(target)) {
    inds <- apply(target, 1, function(r) which.min(sp::spDistsN1(pts = as.matrix(points),
                                                                 pt = c(r[[1]], r[[2]]))))
    return(points[inds, ])
  }
}

yday_366 <- function(date) {
  temp_doy <- lubridate::yday(date)
  temp_leap <- lubridate::leap_year(date)
  temp_doy[(!is.na(temp_doy)) & temp_doy > 59 & (!temp_leap)] <- 1 + temp_doy[(!is.na(temp_doy)) & temp_doy > 59 & (!temp_leap)]
  return(temp_doy)
}

display_daily <- function(Datain,Stations,Variables,option=1,Years,Misscode,Tracecode,Zerocode,Fileout=NA,monstats=c("min","mean","median","max","IQR","sum")){
  #-----------------------------------------------------------------------#
  # Helen Greatrex, SSC, SSC-RCLIM package                                #
  #                                                                       #
  # DisplayDaily: This function displays daily meteorological data in an  #
  #               easy to view format                                     #
  #                                                                       #
  # INPUTS                                                                #
  # Datain   :This is the name of the data.frame containing your data     #
  #           The format of this data.frame is discussed in the details   #
  #           section of this help file                                   #
  #                                                                       #
  # Stations: This is the name of the weather station(s) you wish to view #
  #           e.g. "KUND0002", or c("Paris","London")                     #
  #           It is case sensitive and the stations should be values      # 
  #           in the station column of your data.frame. This is set by    #
  #           default to all stations in the data.frame                   #
  #                                                                       #
  # Variables:This is name of the column of data you wish to view         #
  #           e.g.Rain, TMin, RH.  It is also case sensitive. This is set # 
  #           default to all variables in the data.frame                  #
  #                                                                       #
  # Option:   Set to 1 if you have a column entitled date, which contains #
  #           dates in the format yyyy-mm-dd.  Set to 2 if you have a     #
  #           column entitled Year, one entitled Month and one entitled   #
  #           Day                                                         #
  #                                                                       #
  # Years:    A number or vector of the year(s) you would like to view    #
  #           This is set by default to all years                         #
  #                                                                       #
  # Months:   A number or vector of the month(s) you would like to view   #
  #           This is set by default to all months. This option has been  #
  #            removed as it does not overlap to new years                #
  #                                                                       #
  # Misscode: This is the how you would like to display missing values    #
  #           It is set by default to "-"                                 #
  #                                                                       #
  # Tracecode:This is the how you would like to display trace rainfall    #
  #           values. It is set by default to "tr"                        #
  #                                                                       #
  # Zerocode: This is the how you would like to display zero rainfall     #
  #           values. It is set by default to "--"                        #
  #                                                                       #
  # File:     Would you like to view the results on the screen or print   #
  #           them to a file.  By default this prints to the screen,      #
  #           If you want to print to a file, set this to a character     #
  #           string with the filename (and path if needed)               #
  #                                                                       #
  # DETAILS                                                               #
  # This function prints out daily data in an easy to view format,        #
  # similar to INSTAT's DisplayDaily button.                              # 
  #                                                                       #
  # Your input data.frame should have one column for each variable, one   #
  # row for each date and include the following columns (case sensitive)  #
  #                                                                       #
  # station - containing your station names or identifier                 #
  #                                                                       #
  # EITHER (set option = 1)                                               #
  # date - containing a list of dates in the format yyyy-mm-dd            #
  # OR (set option = 2)                                                   #
  # Year  - containing a list of years in the format yyyy                 #
  # Month - containing a list of months                                   #
  # Day   - containing a list of days of month (from 1-31)                #
  #                                                                       #
  # A column for the variable you wish to study e.g. Rain, TMax, RH, ...  #
  #                                                                       #
  # Note, an experienced R user might wonder why I have not used the      #
  # options in print as much as I might have.  This is because often      #
  # trace rainfall values need to be displayed as text, which messes up   #
  # the print function's missing/zero value options                       #
  #-----------------------------------------------------------------------#
  
  #-----------------------------------------------------------------------#
  # First, check that appropriate inputs exist                            #
  #-----------------------------------------------------------------------#
  if(missing(Datain)==TRUE){stop("Please enter an input data.frame")}
  if(missing(Stations)==TRUE){Stations <- unique(Datain$station)}
  if(missing(Years)==TRUE){Years <- unique(Datain$Year)}
  if(missing(Variables)==TRUE){
    Variables <- names(Datain)
    Variables <- Variables[!(Variables %in% c("Key","station","date","Year","Month","Day","UserDOY"))]
    if(length(Variables) <= 0){
      stop("Please indicate which column name you would like to display")
    }
  }
  flag=0
  Months <- 1:12
  #-----------------------------------------------------------------------#
  # Loop through the stations                                             #
  #-----------------------------------------------------------------------#
  for (s in 1:length(Stations)){
    #-------------------------------------------------------------------#
    # Does the selected station exist?                                  #
    #-------------------------------------------------------------------#
    tmp <- Datain[(as.character(Datain$station) == Stations[s]),]
    
    #-------------------------------------------------------------------#
    # Stop the function if not                                          #
    #-------------------------------------------------------------------#
    if(length(tmp[,1]) <= 0){
      print(paste("The function could not find your chosen station: '",Stations[s],"'",sep=""),quote=FALSE)
      print("Possible choices in the station column of selected data.frame are",quote=FALSE)
      print(unique(Datain$station))
      stop("Choose a station that exists in your data.frame")
    }
    
    #-------------------------------------------------------------------#
    # Loop through the variables                                        #
    #-------------------------------------------------------------------#
    for(v in 1:length(Variables)){
      #--------------------------------------------------------------#
      # Does the selected variable exist?                            #
      #--------------------------------------------------------------#
      if(!(Variables[v] %in% names(tmp))){
        print(paste("There is no column in your data.frame called '", Variables[v],"'",sep=""),quote=FALSE)
        print("Please choose out of",quote=FALSE)
        print(as.character(names(tmp)))
        stop()			
      }
      
      #-------------------------------------------------------------#
      # Now move on to years. First, fill in any missing dates      #
      #-------------------------------------------------------------#
      if(option == 1){
        startday <- as.Date(paste(min(as.numeric(format.Date(tmp$date,"%Y")),na.rm=TRUE),1,1,sep="-"))
        endday <- as.Date(paste(max(as.numeric(format.Date(tmp$date,"%Y")),na.rm=TRUE),12,31,sep="-"))
        tmp$TEMPDate2 <- as.Date(tmp$date)
        tmp$date <- as.Date(tmp$date)
      }else{
        startday <- as.Date(paste(min(tmp$Year,na.rm=TRUE),1,1,sep="-"))
        endday <-   as.Date(paste(max(tmp$Year,na.rm=TRUE),12,31,sep="-"))
        tmp$TEMPDate2 <- as.Date(paste(tmp$Year,tmp$Month,tmp$Day,sep="-"))
        tmp$date <- as.Date(tmp$date)
      }
      
      newdatain <- data.frame(date =  seq(from=startday,to= endday,by="d"))
      new <- merge(newdatain,tmp,by="date",all.x=TRUE,all.y=TRUE)
      # removes NA rows
      new <- new[!is.na(new$TEMPDate2), ]
      loc <- which(names(new)== Variables[v])
      
      # newdates <- seq(startday,endday,1)
      # new <- as.data.frame(matrix(nrow=length(newdates),ncol=dim(tmp)[2]))
      # names(new) <- names(tmp)
      # new$TEMPDate2 <- newdates
      
      # tmpmatch <- match(tmp$TEMPDate2,new$TEMPDate2)
      # new[tmpmatch,(1:(dim(new)[2]-1))] <- tmp[,(1:(dim(new)[2]-1))]
      
      # new <- m
      
      #-------------------------------------------------------------#
      # Rename things like the station name                         #
      #-------------------------------------------------------------#		 
      new$station <- Stations[s]
      
      #-------------------------------------------------------------#
      # Now start to display the data                               #
      # First, make an output template   					  #
      #-------------------------------------------------------------#
      template <- as.data.frame(matrix(nrow=32,ncol=(1+length(Months))))
      names(template) <- c("MONTH",month.abb[Months])
      template[1,] <- "****"
      template[,1] <- c("DAY",1:31)
      
      statstemplate <- as.data.frame(matrix(ncol=(1+length(Months)),nrow=length(monstats)))
      names(statstemplate) <- names(template)
      statstemplate[,1] <- monstats
      
      
      for(y in 1:length(Years)){
        #----------------------------------------------------------#
        # Print out header                                         #
        #----------------------------------------------------------#
        print(" ",quote=FALSE)
        print(paste("STATION  : '",Stations[s],"'",sep=""),quote=FALSE)
        print(paste("VARIABLE : '",Variables[v],"'",sep=""),quote=FALSE)
        print(paste("YEAR     :",Years[y]),quote=FALSE)
        
        #----------------------------------------------------------#           
        # And set up the template/subset the data for that year    #
        #----------------------------------------------------------#
        out <- template
        outstats <- statstemplate
        tmp2 <- new[as.numeric(format.Date(new$TEMPDate2,"%Y")) == Years[y],]	
        
        #----------------------------------------------------------#
        # Now, put into the data.frame for each month..            #
        #----------------------------------------------------------#
        for(m in 1:length(Months)){	
          dat  <- tmp2[as.numeric(format.Date(tmp2$TEMPDate2,"%m")) == Months[m],]	
          
          for(st in 1:length(monstats)){
            if(length(eval(parse(text=paste("suppressWarnings(",monstats[st],"(dat[,loc],na.rm=TRUE))"))))> 1){
              outstats[,(m+1)] <- "NA"
              if(m == 1){
                print("",quote=FALSE)
                print("*********",quote=FALSE)
                print(paste("The '", monstats[st],"' statistic has more than one value in its output, e.g.",sep=""),quote=FALSE)
                print(eval(parse(text=paste("suppressWarnings(",monstats[st],"(dat[,loc],na.rm=TRUE))"))),quote=FALSE)
                print("therefore it has been ignored",quote=FALSE)
                print("",quote=FALSE)
              }
            }else{
              outstats[st,(m+1)] <- ""
              #outstats[st,(m+1)] <- as.character(sprintf("%12.1f",eval(parse(text=paste("suppressWarnings(",monstats[st],"(dat[,loc],na.rm=TRUE))")))))
              #outstats[st,(m+1)] <- gsub("\\s+","",outstats[st,(m+1)]) 
            }
          }
          # If all the data is NA, set the stats to NA (except summary_count_missing)
          if(length(dat[,loc])==length(which(is.na(dat[,loc]==TRUE)))){
            outstats[which(monstats != "summary_count_missing"),(m+1)] <- "NA"
          }
          # If there's no data, set all the stats to NA
          if(length(dat[,loc]) == 0) outstats[,(m+1)] <- "NA"
          #-------------------------------------------------------#
          # Set to 2 d.p.                                         #
          #-------------------------------------------------------#	
          dat[is.na(dat[,loc])==FALSE,loc] <- as.character(sprintf("%12.2f",dat[is.na(dat[,loc])==FALSE,loc]))
          
          #-------------------------------------------------------#
          # Deal with trace rainfall values                       #
          #-------------------------------------------------------#		 	
          tmptrace <- which(dat[,loc] == "0.01")
          if(length(tmptrace) > 0 && !missing(Tracecode)){
            dat[tmptrace,loc] <- Tracecode	
          }
          
          #-------------------------------------------------------#
          # Set to 1 d.p. and remove spaces                       #
          #-------------------------------------------------------#		 	
          dat[,loc] <- substr(dat[,loc],1,nchar(dat[,loc])-1)
          dat[,loc] <- gsub("\\s+","",dat[,loc]) 
          
          #-------------------------------------------------------#
          # Deal with missing values (exc non existant dates)     #
          #-------------------------------------------------------#		 	
          tmpmiss <- which(is.na(dat[,loc]) == TRUE)
          if(length(tmpmiss) > 0 && !missing(Misscode)){
            dat[tmpmiss,loc] <- Misscode	
          }	
          
          #-------------------------------------------------------#
          # Deal with zero values (exc non existant dates)        #
          #-------------------------------------------------------#
          tmpzero <- which(dat[,loc] == "0.0")
          if(length(tmpzero) > 0 && !missing(Zerocode)){
            dat[tmpzero,loc] <- Zerocode	
          }	
          tmpzero <- which(outstats[,(m+1)] == "0.0")
          if(length(tmpzero) > 0){
            outstats[tmpzero,(m+1)] <- "0"	
          }	
          #-------------------------------------------------------#
          # And put into the display data.frame                   #
          #-------------------------------------------------------#
          dat$TEMPday <- as.numeric(format.Date(dat$TEMPDate2,"%d"))
          out[2:32,(m+1)] <- dat[match(1:31,dat$TEMPday),loc]
        }#m
        
        #----------------------------------------------------------#
        # Add the statistics onto the main data.frame              #
        #----------------------------------------------------------#				    
        out[dim(out)[1]+1,] <- "____"
        out[dim(out)[1],1] <- "STATS"
        # temp fix to shorten display name of summary function
        outstats[,1][outstats[,1] == "summary_count_missing"] <- "n_miss"
        out <- rbind(out,outstats)
        
        #----------------------------------------------------------#
        # Set non existant dates to be blank and display           #
        #----------------------------------------------------------#				    
        out[is.na(out) == TRUE] <- " "
        print(out,quote=FALSE,justify="right",row.names=F,digits=1)
        print(" ",quote=FALSE)
        
        #----------------------------------------------------------#
        # Write to file if requested                               #
        #----------------------------------------------------------#
        if(is.na(Fileout) == FALSE){
          if(flag==0){
            fileoutput <- file(Fileout, "w")
            flag = 1
            cat(paste("STATION  : '",Stations[s],"'",sep=""),"\n",file=fileoutput,sep="")
            cat(paste("VARIABLE : '",Variables[v],"'",sep=""),"\n",file=fileoutput,sep="")
            cat(paste("YEAR     :",Years[y]),"\n",file=fileoutput,sep="")
            cat(paste(names(out)),"\n",file=fileoutput,sep="\t")
            write.table(out,file=fileoutput,sep="\t",append=TRUE,quote=FALSE,col.names=FALSE,row.names=FALSE)
          }else{
            cat("\n",paste("STATION  : '",Stations[s],"'",sep=""),"\n",file=fileoutput,sep="")
            cat(paste("VARIABLE : '",Variables[v],"'",sep=""),"\n",file=fileoutput,sep="")
            cat(paste("YEAR     :",Years[y]),"\n",file=fileoutput,sep="")
            cat(paste(names(out)),"\n",file=fileoutput,sep="\t")
            write.table(out,file=fileoutput,sep="\t",append=TRUE,quote=FALSE,col.names=FALSE,row.names=FALSE)
          }
        }
        rm(tmp2)
        rm(out)
      }#y
    }#v
  }#s 
  if(is.na(Fileout) == FALSE){
    close(fileoutput)
  }
} #end of function

#################################################################################

.spells <- function(x, initial_value = NA_real_) {
  y <- mat.or.vec(length(x), 1)
  if(length(x) > 0) {
    y[1] <- dplyr::if_else(x[1] == 0, 0, initial_value + 1)
    if(length(x) > 1) {
      for(i in 2:length(x)) {
        y[i] <- dplyr::if_else(x[i] == 0, 0, y[i-1] + 1)
      }
    }
  }
  return(y)
}

skable <- function(kable_input) {
  kable_input %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                                full_width = FALSE)
}

dekad <- function(date) {
  temp_dekade <- 3 * (lubridate::month(date)) - 2 + (lubridate::mday(date) > 10) + (lubridate::mday(date) > 20)
  return(temp_dekade)
}

pentad <- function(date) {
  temp_pentad <- 6*(lubridate::month(date)) - 5 + (lubridate::mday(date) > 5) + (lubridate::mday(date) > 10) + (lubridate::mday(date) > 15) + (lubridate::mday(date) > 20) + (lubridate::mday(date) > 25)
  return(temp_pentad)	
}