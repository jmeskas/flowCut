#####################################################################
# flowCut: Precise and Accurate Automated Removal of Outlier Events #
#  and Flagging of Files Based on Time Versus Fluorescence Analysis #
# Authors: Justin Meskas and Sherrie Wang                           #
#####################################################################
flowCut <- function(f,
                    Segment=500,
                    Channels=NULL,
                    Directory=NULL,
                    FileID=NULL,
                    Plot="Flagged Only",
                    MaxContin=0.1,
                    MeanOfMeans=0.13,
                    MaxOfMeans=0.15,
                    MaxValleyHgt=0.1,
                    MaxPercCut=0.3,
                    LowDensityRemoval=0.1,
                    GateLineForce=NULL,
                    UseOnlyWorstChannels=FALSE,
                    AmountMeanRangeKeep=1,
                    AmountMeanSDKeep=2,
                    PrintToConsole=FALSE,
                    AllowFlaggedRerun = FALSE,
                    Verbose=FALSE
                    ){

    start0 <- Sys.time()
    resTable <- matrix("", 17, 1)
    rownames(resTable) <-
        c("Is it monotonically increasing in time",
        "Largest continuous jump",
        "Continuous - Pass",
        "Mean of % of range of means divided by range of data",
        "Mean of % - Pass",
        "Max of % of range of means divided by range of data",
        "Max of % - Pass",
        "Has a low density section been removed",
        "% of low density removed",
        "How many segments have been removed",
        "% of events removed from segments removed",
        "Worst channel",
        "% of events removed",
        "FileID",
        "Type of Gating",
        "Was the file run twice",
        "Has the file passed"
        )
    resTable["Was the file run twice",] <- "No"
    # Creating a directory if none is specified
    if (is.null(Directory)){ Directory <- paste0(getwd(), "/flowCut") }
    # Creating a FileID if none is specified
    if (is.null(FileID)) {
        FileID <- Sys.time()
        FileID <- substr(FileID, start = 1, stop = 19)
        FileID <- gsub("-", "_",  FileID)
        FileID <- gsub(":", "_",  FileID)
        FileID <- gsub(" ", "__", FileID)
        # samp <- sample(1:9999, 1)
        # FileID <- paste0(rep(0, 4-nchar(samp)), samp)
        if (Verbose == TRUE){ cat(paste0("The FileID is: ", FileID, "\n"))}
    }
    resTable["FileID", ] <- FileID

    if (class(f) != "flowFrame"){
        message("f must be a flowFrame")
        resTable["Has the file passed", ] <- "f must be a flowFrame"
        return(list(frame=f, ind=NULL, data=resTable, worstChan=NULL))
    }
    if (!is.numeric(Segment) || (Segment <= 0)){
        message("Segment must be a number larger than 0.")
        resTable["Has the file passed", ] <- "Segment must be a number larger than 0."
        return(list(frame=f, ind=NULL, data=resTable, worstChan=NULL))
    }

    if (nrow(f) <= 3*Segment){ # deGate requires count.lim = 3
        message(paste0("Either your Segment size is too large or your number",
                " of cells is too small."))
        resTable["Has the file passed", ] <- "Either your Segment size is too large or your number of cells is too small."
        return(list(frame=f, ind=NULL, data=resTable, worstChan=NULL))
    }
    if (!is.character(Directory)){
        message("Directory must be a character.")
        resTable["Has the file passed", ] <- "Directory must be a character."
        return(list(frame=f, ind=NULL, data=resTable, worstChan=NULL))
    }
    if (!is.character(FileID) && !is.numeric(FileID)){
        message("FileID must be a character or a number.")
        resTable["Has the file passed", ] <- "FileID must be a character or a number."
        return(list(frame=f, ind=NULL, data=resTable, worstChan=NULL))
    }
    if ( !(Plot == "None" || Plot == "All" || Plot == "Flagged Only") ){
        message(paste0("Plot must be a character with one of the following",
                       " options: \'None\', \'All\', or \'Flagged Only\'."))
        resTable["Has the file passed", ] <- "Plot must be a character with one of the following options: \'None\', \'All\', or \'Flagged Only\'."
        return(list(frame=f, ind=NULL, data=resTable, worstChan=NULL))
    }
    if (!is.numeric(MaxContin) || (MaxContin < 0)){
        message("MaxContin must be a number larger than 0.")
        resTable["Has the file passed", ] <- "MaxContin must be a number larger than 0."
        return(list(frame=f, ind=NULL, data=resTable, worstChan=NULL))
    }
    if (!is.numeric(MeanOfMeans) || (MeanOfMeans < 0)){
        message("MeanOfMeans must be a number larger than 0.")
        resTable["Has the file passed", ] <- "MeanOfMeans must be a number larger than 0."
        return(list(frame=f, ind=NULL, data=resTable, worstChan=NULL))
    }
    if (!is.numeric(MaxOfMeans) || (MaxOfMeans < 0)){
        message("MaxOfMeans must be a number larger than 0.")
        resTable["Has the file passed", ] <- "MaxOfMeans must be a number larger than 0."
        return(list(frame=f, ind=NULL, data=resTable, worstChan=NULL))
    }
    if (!is.numeric(MaxValleyHgt) || (MaxValleyHgt < 0) || (MaxValleyHgt > 1)){
        message("MaxValleyHgt must be a number between 0 and 1.")
        resTable["Has the file passed", ] <- "MaxValleyHgt must be a number between 0 and 1."
        return(list(frame=f, ind=NULL, data=resTable, worstChan=NULL))
    }
    if (!is.numeric(MaxPercCut) || (MaxPercCut < 0) || (MaxPercCut > 1)){
        message("MaxPercCut must be a number between 0 and 1.")
        resTable["Has the file passed", ] <- "MaxPercCut must be a number between 0 and 1."
        return(list(frame=f, ind=NULL, data=resTable, worstChan=NULL))
    }
    if (!is.numeric(LowDensityRemoval) || (LowDensityRemoval < 0) ||
        (LowDensityRemoval > 1)){
        message("LowDensityRemoval must be a number between 0 and 1.")
        resTable["Has the file passed", ] <- "LowDensityRemoval must be a number between 0 and 1."
        return(list(frame=f, ind=NULL, data=resTable, worstChan=NULL))
    }
    if (!is.null(GateLineForce) && !is.numeric(GateLineForce) ){
        message("GateLineForce must be numeric.")
        resTable["Has the file passed", ] <- "GateLineForce must be numeric."
        return(list(frame=f, ind=NULL, data=resTable, worstChan=NULL))
    }
    if (!is.logical(UseOnlyWorstChannels)){
        message("UseOnlyWorstChannels must be a logical (boolean).")
        resTable["Has the file passed", ] <- "UseOnlyWorstChannels must be a logical (boolean)."
        return(list(frame=f, ind=NULL, data=resTable, worstChan=NULL))
    }
    if ( AmountMeanRangeKeep %% 1 != 0 ){
        message("AmountMeanRangeKeep must be an integer.")
        resTable["Has the file passed", ] <- "AmountMeanRangeKeep must be an integer."
        return(list(frame=f, ind=NULL, data=resTable, worstChan=NULL))
    }
    if ( AmountMeanSDKeep %% 1 != 0 ){
        message("AmountMeanSDKeep must be an integer.")
        resTable["Has the file passed", ] <- "AmountMeanSDKeep must be an integer."
        return(list(frame=f, ind=NULL, data=resTable, worstChan=NULL))
    }
    if (!is.logical(PrintToConsole)){
        message("PrintToConsole must be a logical (boolean).)")
        resTable["Has the file passed", ] <- "PrintToConsole must be a logical (boolean).)"
        return(list(frame=f, ind=NULL, data=resTable, worstChan=NULL))
    }
    if (!is.logical(Verbose)){
        message("Verbose must be a logical (boolean).")
        resTable["Has the file passed", ] <- "Verbose must be a logical (boolean)."
        return(list(frame=f, ind=NULL, data=resTable, worstChan=NULL))
    }

    #### Extracting the location of channels where flowCut will be applied ####
    t.name <- f@parameters@data$name
    t.desc <- f@parameters@data$desc
    FSC.loc  <- sort(unique(c(
                    grep("fsc", tolower(t.name)),
                    grep("fs lin", tolower(t.name)),
                    grep("*FS", t.name)
                )))
    names( FSC.loc) <- NULL
    SSC.loc  <- sort(unique(c(grep("ssc", tolower(t.name)),
                    grep("ss lin", tolower(t.name)),
                    grep("*SS", t.name)
                )))
    names( SSC.loc) <- NULL

    Time.loc <- sort(unique(c(grep("time",  tolower(t.name)),
                              grep("time",  tolower(t.desc)),
                              grep("hdr-t", tolower(t.name))
                )))

    all.Time.loc <- Time.loc
    # if (length(Time.loc) >= 2){
    #   Time.loc <- Time.loc[1] # default to take the first.
    # }
    # names(Time.loc) <- NULL
    #

    if (length(all.Time.loc) >= 2){ # Multiple time channels
        message(paste0("This file has ", length(all.Time.loc),
                       " time channels. flowCut has selected to use ",
                       t.name[Time.loc], " - ", t.desc[Time.loc], "."))
        Time.loc <- Time.loc[1] # default to take the first.

        f@parameters@data$name[all.Time.loc[which(all.Time.loc != Time.loc)]] <-
            paste0(f@parameters@data$name[all.Time.loc[which(all.Time.loc != Time.loc)]], "-Removed")
        colnames(f@exprs)[all.Time.loc[which(all.Time.loc != Time.loc)]] <-
            paste0(colnames(f@exprs)[all.Time.loc[which(all.Time.loc != Time.loc)]], "-Removed")
    }

    Extra.loc <-  c(
        grep( "pulse",           tolower(t.name)),
        grep( "width",           tolower(t.name)),
        grep( "length",          tolower(t.name)),
        grep( "count",           tolower(t.name)),
        grep( "sort classifier", tolower(t.name)),
        grep( "event",           tolower(t.name)),
        grep( "phenograph",      tolower(t.name)),
        grep( "barcode",         tolower(t.name))
    )
    names(Extra.loc) <- NULL
    Extra.loc <- unique(Extra.loc)
    if (length(Extra.loc) >= 1 && Verbose == T ){
        cat(paste0("Channels ", paste0(Extra.loc, collapse = ", "), " are removed as they are not channels that need to be analyzed.\n"))
    }

    NoVariation <- NULL
    for (NoVar in 1:length(colnames(f))){
        if (sd(f@exprs[,NoVar], na.rm = T) == 0){
            NoVariation <- c(NoVariation, NoVar)
        }
    }
    names(NoVariation) <- NULL
    if (length(NoVariation) >= 1 && Verbose == T ){
        message(paste0("Channels ", paste0(NoVariation, collapse = ", "), " have no variation and have been removed from the analysis."))
    }



    MonotonicWithTime <- NULL
    for (MonoChan in 1:length(colnames(f))){
        if (all(f@exprs[ , NoVar] == cummax(f@exprs[ , NoVar])) == TRUE){
            MonotonicWithTime <- c(MonotonicWithTime, NoVar)
        }
    }
    names(MonotonicWithTime) <- NULL
    MonotonicWithTime <- sort(unique(MonotonicWithTime))

    if ( length(which(match(all.Time.loc, MonotonicWithTime, nomatch=0) >= 1)) >= 1){
        MonotonicWithTime <- MonotonicWithTime[-match(all.Time.loc, MonotonicWithTime, nomatch=0)]
    }

    if (length(MonotonicWithTime) >= 1 && Verbose == T ){
        message(paste0("Channels ", paste0(MonotonicWithTime, collapse = ", "), " are monotonically increasing in time and have been removed from the analysis."))
    }




    if(length(which(NoVariation == Time.loc)) >= 1){
        message("Your time channel has no variation.")
        f@parameters@data$name[Time.loc] <- paste0(f@parameters@data$name[Time.loc], "-Removed")
        colnames(f@exprs)[Time.loc] <- paste0(colnames(f@exprs)[Time.loc], "-Removed")
        Time.loc <- NULL

        # if(length(which(NoVariation != all.Time.loc)) >= 1){
        if(length(which(is.na(match(all.Time.loc,NoVariation)))) >= 1){ # does something exist in all.Time.loc that isn't in NoVariation

            message("The first time channel will be replaced by the second time channel")
            Time.loc <- all.Time.loc[-which(all.Time.loc == NoVariation)][1]
            f@parameters@data$name[Time.loc] <- "Time"
            f@parameters@data$name[all.Time.loc[which(all.Time.loc != Time.loc)]] <-
                paste0(f@parameters@data$name[all.Time.loc[which(all.Time.loc != Time.loc)]], "-Removed")
            colnames(f@exprs)[all.Time.loc[which(all.Time.loc != Time.loc)]] <-
                paste0(colnames(f@exprs)[all.Time.loc[which(all.Time.loc != Time.loc)]], "-Removed")
        }
    }

    #### Creating a time channel if none is specified #########################
    if (length(Time.loc) == 0){
        message("Your data does not have a time channel. flowCut will",
                " create one, but now flowCut will not be as fully",
                " functional as it could be. Consider recording the time",
                " for future projects.")
        f@exprs <- cbind(f@exprs, 1:nrow(f))
        colnames(f@exprs)[length(colnames(f))+1] <- "Time"
        f@parameters@data <- rbind(
            f@parameters@data,
            c("Time", "Time", 262144, -111, 262143)
        )
        rownames(f@parameters@data)[length(colnames(f))] <-
            paste0("$P", length(colnames(f)))
        f@description[paste0("P", length(colnames(f)), "DISPLAY")] <- "LIN" # "LOG"
        f@description[paste0("flowCore_$P", length(colnames(f)), "Rmax")] <- 262143
        f@description[paste0("flowCore_$P", length(colnames(f)), "Rmin")] <- 0
        f@description[paste0("P", length(colnames(f)), "B")] <- "0"
        f@description[paste0("P", length(colnames(f)), "R")] <- "262143"
        f@description[paste0("P", length(colnames(f)), "N")] <- "Time"
        f@description[paste0("P", length(colnames(f)), "G")] <- "1"
        f@description[paste0("P", length(colnames(f)), "E")] <- "0,0"


        Time.loc <- length(colnames(f))
    }

    if (length(c(FSC.loc, SSC.loc)) == 0 ){
        message("No FCS or SSC channels found.")
        # return(list(frame=f, ind=NULL, data=resTable, worstChan=NULL))
    }

    range_of_time <- max(f@exprs[,Time.loc])-min(f@exprs[,Time.loc])

    Time_test_passes <- TRUE

    if ( range_of_time != 0){

        # is the mean and midpoint similar?
        uniformity_in_time_test <- abs(mean(f@exprs[,Time.loc]) - (range_of_time/2+min(f@exprs[,Time.loc]))) / range_of_time

        # print(uniformity_in_time_test)

        # dividing_points <- seq(min(f@exprs[,Time.loc]), max(f@exprs[,Time.loc]), length.out = 11)
        # uniformity_in_time_test_2 <- sapply (1:10, function(x) {
        #     length(intersect(
        #         which(f@exprs[,Time.loc] <= dividing_points[x+1]),
        #         which(f@exprs[,Time.loc] >= dividing_points[x]))
        #     )
        # })

        # This number needs to be set better. Used to be 0.2, moved to 0.22 to work with a particular project.
        if ( uniformity_in_time_test >= 0.22 ){
            message("The time channel does not appear to be distributed like an expected time channel would be.")
            Time_test_passes <- FALSE
        }

        # if ( (min(uniformity_in_time_test_2) <= 0.05*max(uniformity_in_time_test_2)) ){
            # message("The 2 time channel does not appear to be distributed like an expected time channel would be.")
            # Time_test_passes <- FALSE
        # }

        # is there a bunch of repeats?
        uniformity_in_time_test_3 <- max(table(f@exprs[,Time.loc]))

        # print(uniformity_in_time_test_3)

        if ( uniformity_in_time_test_3 >= 0.05*nrow(f) ){
            message("There appears to be an overwhelming amount of repeated time values.")
            Time_test_passes <- FALSE
        }
    } else { # if all the same time value, messes up other time tests
        message("Range of the time channel is 0.")
        Time_test_passes <- FALSE
    }

    if(Time_test_passes == FALSE){
        message("Time test(s) failed.")
        resTable["Has the file passed", ] <- "Time test(s) failed."
        return(list(frame=f, ind=NULL, data=resTable, worstChan=NULL))
    }

    # channels to clean
    CleanChan.loc <- (1:ncol(f))[-c(FSC.loc, SSC.loc, Time.loc, all.Time.loc, Extra.loc, NoVariation, MonotonicWithTime)]

    if (length(CleanChan.loc) == 0 ){
        message("No marker channels to run flowCut on.")
        resTable["Has the file passed", ] <- "No marker channels to run flowCut on."
        return(list(frame=f, ind=NULL, data=resTable, worstChan=NULL))
    }

    # Use all non scatter/time channels unless otherwise specified
    if (!is.null(Channels)){
        if (all(is.character(Channels))){
            Channels <- sort(unique(sapply(1:length(Channels), function(x) {
                grep(tolower(Channels[x]), tolower(f@parameters@data$desc))
            })))
        }

        # Forces the user not to use FSC, SSC and Time
        CleanChan.loc <- intersect(CleanChan.loc, Channels)
    }

    ind.removed <- NA
    f.org <- f

    #### Test if the file is monotonic ########################################
    if ( all(f@exprs[ , Time.loc] == cummax(f@exprs[ , Time.loc])) == FALSE ){
        message("The flow frame is not monotonically increasing in time")
        resTable["Is it monotonically increasing in time", ] <- "F"
    } else {
        resTable["Is it monotonically increasing in time", ] <- "T"
    }

    #### Remove low density sections ##########################################
    res.temp <- removeLowDensSections(f, Segment, LowDensityRemoval, Verbose=Verbose, Time.loc=Time.loc)
    f <- res.temp$frame; removeIndLowDens <- res.temp$rem.ind; remove(res.temp)
    ifelse(length(removeIndLowDens) >= 1 ,
        resTable["Has a low density section been removed", ] <- "T",
        resTable["Has a low density section been removed", ] <- "F")
    resTable["% of low density removed", ] <-
        as.character( round(length(removeIndLowDens) /
            nrow(f.org), digits=4) * 100 )

    #### 1st Time: Calculate means and which segments will be removed #########
    res.temp <- calcMeansAndSegmentsRemoved(
        f=f,
        Segment=Segment,
        CleanChan.loc=CleanChan.loc,
        FirstOrSecond="First",
        MaxValleyHgt=MaxValleyHgt,
        MaxPercCut=MaxPercCut,
        MaxContin=MaxContin,
        MeanOfMeans=MeanOfMeans,
        MaxOfMeans=MaxOfMeans,
        GateLineForce=GateLineForce,
        UseOnlyWorstChannels=UseOnlyWorstChannels,
        AmountMeanRangeKeep=AmountMeanRangeKeep,
        AmountMeanSDKeep=AmountMeanSDKeep,
        Verbose=Verbose,
        Time.loc=Time.loc
    )

    deletedSegments1 <- res.temp$deletedSegments;
    quantiles1       <- res.temp$quantiles;
    storeMeans1      <- res.temp$storeMeans
    meanRangePerc1   <- res.temp$meanRangePerc;
    timeCentres1     <- res.temp$timeCentres;
    typeOfGating     <- res.temp$typeOfGating
    densityGateLine  <- res.temp$densityGateLine;
    cellDelete2      <- res.temp$cellDelete2;
    choosenChans     <- res.temp$choosenChans
    remove(res.temp)

    #### Now delete segments that are statistically different #################
    removed.ind <- NULL
    totalNumSeg <- floor(nrow(f@exprs)/Segment)

    if ( length(deletedSegments1) == 0 ){ # If nothing is deleted
        if (Verbose == TRUE){
            cat("None deleted from flowCut segment removal.\n")
        }
        resTable["How many segments have been removed", ] <- as.character(0)
    } else {
        deletedSegments1 <-
            sort( unique( deletedSegments1 ), decreasing = FALSE )

        # Turn segments removed into ranges for printing to console
        del.seg.list <-
            split(deletedSegments1, cumsum(c(1, diff(deletedSegments1) != 1)))
        print.segs.rem <- sapply(1:length(del.seg.list), function(x) {
            if( length(del.seg.list[[x]]) >= 2) {
                paste0(del.seg.list[[x]][1], "-",
                    del.seg.list[[x]][length(del.seg.list[[x]])]
                )
            } else {
                paste0(del.seg.list[[x]][1])
            }
        } )

        if (Verbose == TRUE){
            cat(paste0("Removing segments ",
                paste0(print.segs.rem, collapse = ", "), " out of ",
                    totalNumSeg, " segments.\n"
                )
            )
        }
        resTable["How many segments have been removed", ] <-
            as.character( length(deletedSegments1) )
        for (n in 1:length(deletedSegments1)){
            if (deletedSegments1[n] == totalNumSeg)
                removed.ind <-
                    c(
                        removed.ind,
                        (Segment*(deletedSegments1[n]-1)+1):nrow(f@exprs)
                    )
            if (deletedSegments1[n] != totalNumSeg)
                removed.ind <-
                    c(removed.ind, Segment*(deletedSegments1[n]-1)+(1:Segment))
        }
        f@exprs <- f@exprs[-removed.ind, ]
    }
    resTable["% of events removed from segments removed", ] <-
        as.character( round(length(removed.ind) / nrow(f.org), digits=4) * 100)

    #### 2nd Time: Calculate means and which segments will be removed #########
    # This time we are only interested in calculations of quantiles, means,
    #  mean range percentage and we do not want to do any deletion
    res.temp <- calcMeansAndSegmentsRemoved(
        f=f,
        Segment=Segment,
        CleanChan.loc=CleanChan.loc,
        FirstOrSecond="Second",
        MaxValleyHgt=MaxValleyHgt,
        MaxPercCut=MaxPercCut,
        MaxContin=MaxContin,
        MeanOfMeans=MeanOfMeans,
        MaxOfMeans=MaxOfMeans,
        GateLineForce=GateLineForce,
        UseOnlyWorstChannels=UseOnlyWorstChannels,
        AmountMeanRangeKeep=AmountMeanRangeKeep,
        AmountMeanSDKeep=AmountMeanSDKeep,
        Verbose=Verbose,
        Time.loc=Time.loc
    )
    quantiles2 <- res.temp$quantiles;
    storeMeans2 <- res.temp$storeMeans;
    meanRangePerc2 <- res.temp$meanRangePerc;
    timeCentres2 <- res.temp$timeCentres;
    remove(res.temp)

    #### Continuous Check #####################################################
    # Check if there are sudden changes in the mean for the neighbouring
    #  segments. If the difference of the mean of a given segment to the
    #  mean of an adjacent segment divided by the difference of the 98th
    #  percentile and the 2nd percentile for the whole file is larger than
    #  a critical value, then the file is flagged.

    maxDistJumped <- rep(0, max(CleanChan.loc))
    for ( j in CleanChan.loc ){
        temp.vect <- rep(0, length(storeMeans2[[j]])-1)
        for ( l in 1:(length(storeMeans2[[j]])-1) ){
            temp.vect[l] <- abs(storeMeans2[[j]][l]-storeMeans2[[j]][[l+1]]) /
                (quantiles1[[j]]["98%"]-quantiles1[[j]]["2%"])
        }
        maxDistJumped[j] <- max(temp.vect)
    }
    resTable["Largest continuous jump", ] <-
        as.character( round(max(maxDistJumped, na.rm = TRUE), digits=3) )
    if ( resTable["Largest continuous jump", ] >= as.numeric(MaxContin)) {
        resTable["Continuous - Pass", ] <- "F"
        if ( Verbose == TRUE){
            message("The file has been flagged. The largest continuous jump ",
                "was larger than ", MaxContin*100, "% of the range of the ",
                "2-98 percentile of the full data.")
        }
    } else {
        resTable["Continuous - Pass", ] <- "T"
    }

    #### Mean of Range of Means Divided by Range of Data ######################
    # If there is a large gradual change of the means of the fluorescence in
    #  all channels, then the file is flagged.
    resTable["Mean of % of range of means divided by range of data", ] <-
        as.character( round(mean(meanRangePerc2, na.rm = TRUE), digits=3) )
    if ( resTable["Mean of % of range of means divided by range of data", ] >=
        as.numeric(MeanOfMeans)) {
        if ( Verbose == TRUE ){
            message("The file has been flagged. The means differ more than ",
                MeanOfMeans*100, "% of the range of the 2-98 percentile of ",
                "the full data.")
        }
        resTable["Mean of % - Pass", ] <- "F"
    } else {
        resTable["Mean of % - Pass", ] <- "T"
    }

    #### Max of Range of Means Divided by Range of Data ####
    # If there is a very large gradual change of the means of the fluorescence
    #  in at least one channels, then the file is flagged.
    resTable["Max of % of range of means divided by range of data", ] <-
        round(max(meanRangePerc2, na.rm = TRUE), digits=3)
    # use 1 because we want the worst marker before any corrections.
    worstChan <- min(which(meanRangePerc1 == max(meanRangePerc1, na.rm=TRUE)))
    names.worschan <- f@parameters@data$name[worstChan];
    names(names.worschan) <- NULL
    resTable["Worst channel", ] <- names.worschan
    if ( resTable["Max of % of range of means divided by range of data", ] >=
        MaxOfMeans) {
        if ( Verbose == TRUE ){
            message("The file has been flagged. The max ranged means differ ",
            "more than ", MaxOfMeans*100, "% of the range of the 2-98 ",
            "percentile of the full data.")
        }
        resTable["Max of % - Pass", ] <- "F"
    } else {
        resTable["Max of % - Pass", ] <- "T"
    }

    #### Organize the indices that have been removed ##########################
    if (  is.null(removed.ind ) &&  is.null(removeIndLowDens))
        to.be.removed <- NULL
    if (  is.null(removed.ind ) && !is.null(removeIndLowDens))
        to.be.removed <- removeIndLowDens
    if ( !is.null(removed.ind ) &&  is.null(removeIndLowDens))
        to.be.removed <- removed.ind
    if ( !is.null(removed.ind ) && !is.null(removeIndLowDens)){
        # lowDens was removed first
        temp <- setdiff(1:nrow(f.org), removeIndLowDens)
        to.be.kept <- temp[setdiff(1:length(temp), removed.ind)]
        to.be.removed <- setdiff(1:nrow(f.org), to.be.kept)
    }

    resTable["% of events removed", ] <-
        as.character(round(length(to.be.removed) / nrow(f.org), digits=4)*100)
    if ( "TTTT" == paste0(resTable["Is it monotonically increasing in time", ],
                        resTable["Continuous - Pass", ],
                        resTable["Mean of % - Pass", ],
                        resTable["Max of % - Pass", ]) ){
        resTable["Has the file passed", ] <- "T"
    } else {
        resTable["Has the file passed", ] <- "F"
    }

    # Keep track of the worse channels by using asterisks
    asterisks <- rep("", max(CleanChan.loc))
    if(UseOnlyWorstChannels == TRUE){
        asterisks <- sapply(1:max(CleanChan.loc), function(x) {
            if(length(which(x == choosenChans)) >= 1){
                asterisks <- " *"
            } else {
                asterisks <- ""
            }
            return(asterisks)
        })
    }

    if (Verbose == TRUE){ cat("Type of Gating: ", typeOfGating, ".\n", sep="")}
    resTable["Type of Gating",] <- typeOfGating

    #### Plotting #############################################################
    if (resTable["Is it monotonically increasing in time", ] == "T"){
        PassedMono   <- "T"     } else { PassedMono   <- "F" }
    if (resTable["Continuous - Pass", ]                      == "T"){
        PassedCont   <- "T"     } else { PassedCont   <- "F" }
    if (resTable["Mean of % - Pass", ]                       == "T"){
        PassedMean   <- "T"     } else { PassedMean   <- "F" }
    if (resTable["Max of % - Pass", ]                        == "T"){
        PassedMax    <- "T"     } else { PassedMax    <- "F" }
    if (resTable["Has the file passed", ]                    == "T"){
        FlaggedOrNot <- "Passed"} else { FlaggedOrNot <- "Flagged" }

    # the pngName will be passed through the variable "AllowFlaggedRerun" for the second run to make the second figure had a suffix.
    pngName <- paste0(Directory, "/", gsub(".fcs","",FileID), "_", FlaggedOrNot, "_", PassedMono, PassedCont, PassedMean, PassedMax, ".png")

    if(Plot == "All" || (Plot == "Flagged Only" && FlaggedOrNot == "Flagged")){
        # z1 and z2 are the dimensions of the figure
        z1 <- ceiling(sqrt(length(CleanChan.loc)+2))
        if ( (z1^2-z1) >= (length(CleanChan.loc)+2)){
            z2 <- z1-1
        } else {
            z2 <- z1
        }

        suppressWarnings ( dir.create ( paste0(Directory), recursive = TRUE) )

        if ( AllowFlaggedRerun != T && AllowFlaggedRerun != F && file.exists(AllowFlaggedRerun) ){
            pngName <- gsub(".png", "_2nd_run.png", pngName)
        }

        if ( PrintToConsole == FALSE) {
            CairoPNG ( filename = pngName, width = (z1)*600, height = z2*600)
            par(mfrow=c(z2,z1), mar=c(7,7,4,2), mgp=c(4,1.5,0), oma=c(0,0,5,0))
            cex.size <- 3
        } else {
            par(mfrow=c(z2,z1), mar=c(5,5,4,2), mgp=c(3,1,0), oma=c(0,0,5,0))
            cex.size <- 1.5
        }

        for (x in CleanChan.loc){
            plotDens(f.org,
                    c(Time.loc,x),
                    cex.main=cex.size,
                    cex.lab=cex.size,
                    cex.axis=cex.size,
                    main=paste0(
                        round(meanRangePerc1[x], digits=3), " / ",
                        round(meanRangePerc2[x], digits=3), " (",
                        round(max(maxDistJumped[x]), digits=3), ")",
                        asterisks[x])
                    )

            if ( length(to.be.removed) != 0 )
                points(exprs(f.org)[to.be.removed, c(Time.loc,x)],
                    col=1, pch=".")
            if ( (length(removeIndLowDens) != 0 ) )
                points(exprs(f.org)[removeIndLowDens, c(Time.loc,x)],
                    pch=".", cex=1, col="grey")
            lines(x=(timeCentres1), y=storeMeans1[[x]], cex.main=cex.size,
                cex.lab=cex.size, cex.axis=cex.size, lwd=4, col="deeppink2")
            lines(x=(timeCentres2), y=storeMeans2[[x]], cex.main=cex.size,
                cex.lab=cex.size, cex.axis=cex.size, lwd=4, col="brown")
            abline(h=c(quantiles1[[x]]["98%"], quantiles1[[x]]["2%"]),
                lwd=4, col="chocolate2")
            abline(h=c(quantiles2[[x]]["98%"], quantiles2[[x]]["2%"]),
                lwd=4, col="chocolate4")
        }
        x <- worstChan

        plotDens(f.org,
            c(Time.loc,x),
            cex.main=cex.size,
            cex.lab=cex.size,
            cex.axis=cex.size,
            main=paste0("Worst Channel without indices removed")
        )

        temp <- density(cellDelete2, adjust = 1)
        graphics::plot(temp, cex.main=cex.size, cex.lab=cex.size,
            cex.axis=cex.size, main='Density of summed measures')
        abline(v=densityGateLine, lwd=2)
        title(main=paste0(FileID, " ", FlaggedOrNot, " ", PassedMono,
            PassedCont, PassedMean, PassedMax),
            outer=TRUE, cex.main=cex.size+1
        )
        if ( PrintToConsole == FALSE) {
            dev.off()
        } else {
            par(mfrow=c(1,1), mar=c(5,5,4,2), mgp=c(3,1,0)) # back to default
        }
    }

    if (Verbose == TRUE){
        if(resTable["Has the file passed", ] == "T") {
            cat("File Passed\n")
        } else {
            cat(paste0("The file has been flagged ",
                PassedMono, PassedCont, PassedMean, PassedMax, "\n"))
        }
    }

    if (Verbose == TRUE){
        cat("Cleaning completed in: ", TimePrint(start0), "\n", sep="")
    }

    if (AllowFlaggedRerun == TRUE && resTable["Has the file passed", ] == "F"){
        if (Verbose == TRUE){ cat("Running flowCut a second time.\n")}
        res_flowCut <- flowCut(
            f=f,
            Segment=Segment,
            Channels=Channels,
            Directory=Directory,
            FileID=FileID,
            Plot=Plot,
            MaxContin=MaxContin,
            MeanOfMeans=MeanOfMeans,
            MaxOfMeans=MaxOfMeans,
            MaxValleyHgt=MaxValleyHgt,
            MaxPercCut=MaxPercCut,
            LowDensityRemoval=LowDensityRemoval,
            GateLineForce=GateLineForce,
            UseOnlyWorstChannels=UseOnlyWorstChannels,
            AmountMeanSDKeep=AmountMeanSDKeep,
            AmountMeanRangeKeep=AmountMeanRangeKeep,
            PrintToConsole=PrintToConsole,
            AllowFlaggedRerun=pngName,
            Verbose=Verbose
        )

        indOfInd <- setdiff(1:nrow(f.org), to.be.removed)
        indOfInd <- sort(c(indOfInd[res_flowCut$ind], to.be.removed))

        resTableOfResTable <- res_flowCut$data
        if ( as.numeric(res_flowCut$data["Largest continuous jump", ])
            <  as.numeric(resTable["Largest continuous jump", ])){
            resTableOfResTable["Largest continuous jump", ] <-
                resTable["Largest continuous jump", ]
        }
        if ( as.numeric(res_flowCut$data["Mean of % of range of means divided by range of data", ])
                <  as.numeric(resTable["Mean of % of range of means divided by range of data", ])){
            resTableOfResTable["Mean of % of range of means divided by range of data", ] <-
                resTable["Mean of % of range of means divided by range of data", ]
        }
        if ( as.numeric(res_flowCut$data["Max of % of range of means divided by range of data", ])
                <  as.numeric(resTable["Max of % of range of means divided by range of data", ])){
            resTableOfResTable["Max of % of range of means divided by range of data", ] <-
                resTable["Max of % of range of means divided by range of data", ]
        }
        resTableOfResTable["% of low density removed", ] <-
            as.character( round((nrow(f)*as.numeric(res_flowCut$data["% of low density removed", ])
                + nrow(f.org)*as.numeric(resTable["% of low density removed", ]))/ nrow(f.org), digits=4) )
        resTableOfResTable["How many segments have been removed", ] <-
            as.character(as.numeric(res_flowCut$data["How many segments have been removed", ])
                + as.numeric(resTable["How many segments have been removed", ]))
        resTableOfResTable["% of events removed from segments removed", ] <-
            as.character( round((nrow(f)*as.numeric(res_flowCut$data["% of events removed from segments removed", ])
                + nrow(f.org)*as.numeric(resTable["% of events removed from segments removed", ]))/ nrow(f.org), digits=4) )
        resTableOfResTable["% of events removed", ] <-
            as.character( round((nrow(f)*as.numeric(res_flowCut$data["% of events removed", ])
                + nrow(f.org)*as.numeric(resTable["% of events removed", ]))/ nrow(f.org), digits=4) )
        resTableOfResTable["Was the file run twice",] <- "Yes"
        return(list(frame=res_flowCut$frame, ind=indOfInd, data=resTableOfResTable, worstChan=res_flowCut$worstChan))
    }
    return(list(frame=f, ind=to.be.removed, data=resTable, worstChan=worstChan))
}



################################################
# Remove sections that have a very low density #
################################################
# The sections that have a density of less than LowDensityRemoval
#  (defaulted at 10%) of the maximum density are removed.
# An additional 2.5% of the range of these low density regions on either side
#  is also removed because it is very common to have a burst of events before
#  or after these low density sections.
# Because the density functions indices do not match with the flowFrames
#  indices, density function indices need to convert to flowFrame indices.
removeLowDensSections <- function(f, Segment=500, LowDensityRemoval=0.1, Verbose=FALSE, Time.loc){

    if(LowDensityRemoval == 0){ # dont want to remove low density on the rerun
        return(list(frame=f, rem.ind=NULL))
    }

    # Time.loc <- which(tolower(f@parameters@data$name) == "time")
    # names(Time.loc) <- NULL

    # only plays a role if users are using removeLowDensSections independently
    #  of flowCut
    # if(length(Time.loc) == 0){
    #     message("There is no time channel. removeLowDensSections only works ",
    #             "with a time channel.")
    #     return(list(frame=f, rem.ind=NULL))
    # }

    minTime <- min(f@exprs[ , Time.loc])
    maxTime <- max(f@exprs[ , Time.loc])

    # Change flow frame data into a density
    dens.f <- density(f@exprs[ , Time.loc], n= nrow(f), adjust = 0.1)
    # Extracting indices that has density values below 10% or a user-specified threshold
    low.dens.removeIndLowDens <-
        which(dens.f$y <= LowDensityRemoval*max(dens.f$y))

    # Split the consectutive sections into separate elements in the list
    range.low.dens <- split(
        low.dens.removeIndLowDens,
        cumsum(c(1, diff(low.dens.removeIndLowDens) != 1))
    )

    # If there are indices left in range.low.dens to be removed
    if (length(range.low.dens) != 0 ){

        # change groups of indices to ranges
        range.low.dens <- lapply(1:length(range.low.dens), function(x){
            range(range.low.dens[[x]])
        })
        # Add 2.5% each way to the range.
        range.low.dens <- lapply(1:(length(range.low.dens)), function(x){
            range.temp <- range.low.dens[[x]][2]- range.low.dens[[x]][1]

            # if (range.temp <= 0.25 *(maxTime-minTime)){ # only add buffer if range is less than 25% of the total time range
            #     new.range <- c( max(round(range.low.dens[[x]][1]-0.025*range.temp), 1),
            #                     min(round(range.low.dens[[x]][2]+0.025*range.temp),
            #                         length(dens.f$y))
            #                  )
            # } else {
                new.range <- c( max(round(range.low.dens[[x]][1]), 1),
                                min(round(range.low.dens[[x]][2]),
                                    length(dens.f$y))
                             )
            # }
            return(new.range)
        })

        # Change density function indices to time coordinates
        range.low.dens <- lapply(1:length(range.low.dens), function(x){
            c(dens.f$x[range.low.dens[[x]][1]],
              dens.f$x[range.low.dens[[x]][2]]
            )
        })

        removeIndLowDens <- NULL

        # Change time coordinates to flowframe indices, removeIndLowDens
        #  contains flowframe indices to be removed
        for ( b2 in 1:length(range.low.dens) ){
            removeIndLowDens <- c(
                removeIndLowDens,
                intersect(
                    which(f@exprs[ , Time.loc] >= range.low.dens[[b2]][1]),
                    which(f@exprs[ , Time.loc] <= range.low.dens[[b2]][2])
                )
            )
        }
        removeIndLowDens <- unique(removeIndLowDens)

        if (length(removeIndLowDens) == 0 ){
            if (Verbose == TRUE){
                cat("None deleted from flowCut low dens removal.\n")
            }
        } else if (length(removeIndLowDens) == nrow(f)){
            if (Verbose == TRUE){
                cat("Low density removal removed all events. Probably a spike in the time channel. Therefore, removing no events.\n")
            }
            removeIndLowDens <- NULL
        } else if (length(removeIndLowDens) > 0.5*nrow(f)){
            if (Verbose == TRUE){
                cat("Low density removal removed more than half of the events. Probably a spike in the marker channels. Therefore, removing no events.\n")
            }
            removeIndLowDens <- NULL
        } else if ( (nrow(f) - length(removeIndLowDens)) < Segment*3 ){
            if (Verbose == TRUE){
                cat("Low density removal removed too many events. If we allowed removal of all low density events, then there",
                    "would be less than Segment*3 events and then many errors would occur in flowCut. Therefore, we are not removing any.\n")
            }
            removeIndLowDens <- NULL
        }  else {
            temp.text <- range.low.dens # Create temp.text for printing to console purposes

            # Check if the ranges of time coordinates exceed the maxTime and minTime, if so, update the end points
            if( temp.text[[length(temp.text)]][1] < maxTime &&
                temp.text[[length(temp.text)]][2] > maxTime ){
                temp.text[[length(temp.text)]][2] <- maxTime
            }
            if( temp.text[[length(temp.text)]][1] > maxTime &&
                temp.text[[length(temp.text)]][2] > maxTime ){
                temp.text[[length(temp.text)]] <- NULL
            }

            if ( temp.text[[1]][1] < minTime && temp.text[[1]][2] > minTime ){
                temp.text[[1]][1] <- minTime
            }
            if ( temp.text[[1]][1] < minTime && temp.text[[1]][2] < minTime ){
                temp.text[[1]] <- NULL
            }

            for( p1 in 1:length(temp.text)){
                temp.text[[p1]] <-
                    paste(round(temp.text[[p1]], digits = 2), collapse=" to ")
            }
            # Converts to text-form for output
            temp.text <- paste(temp.text, collapse = ", ")

            if (Verbose == TRUE){
                cat(paste0("Removing low density ranges ", temp.text, ".\n"))
            }
            f@exprs <- f@exprs[-removeIndLowDens, ]
        }
    } else { # If there are no indices to be removed
        if (Verbose == TRUE){
            cat("None deleted from flowCut low dens removal.\n")
        }
        removeIndLowDens <- NULL
    }
    return(list(frame=f, rem.ind=removeIndLowDens))
}



######################################################
# Calculate means and which segments will be removed #
######################################################
# All events are broken down into segments for analysis
# Eight measures of each segment (mean, median, 5th, 20th, 80th and 95th
#  percentile, second moment(variation) and third moment (skewness)) are
#  calculated. The regions where these eight measures are significantly
#  different from the rest will be removed.
calcMeansAndSegmentsRemoved <- function(
    f,
    Segment,
    CleanChan.loc,
    FirstOrSecond,
    MaxValleyHgt,
    MaxPercCut,
    MaxContin,
    MeanOfMeans,
    MaxOfMeans,
    GateLineForce,
    UseOnlyWorstChannels,
    AmountMeanRangeKeep,
    AmountMeanSDKeep,
    Verbose,
    Time.loc
    ){

    # f is the flow frame
    # FirstOrSecond - parameter used to determine what calculation or process
    #  will be performed

    deletedSegments <- NULL
    meanRangePerc   <- NULL
    timeCentres     <- NULL
    cellDelete <- list()
    storeMeans <- list()
    quantiles  <- list()
    totalNumSeg <- floor(nrow(f@exprs)/Segment)

    maxDistJumped <- rep(0, max(CleanChan.loc))
    # Calcuate all means except for the last segment.
    # Each segment containing 500 (Segments = 500) events


    for ( k in 1:(totalNumSeg-1)){
        temp <- f@exprs[Segment*(k-1)+(1:Segment), Time.loc]
        # timeCentres contains the means of each segment
        timeCentres[k] <- mean(temp)
    }
    # Calculating the mean of the last segment
    k <- totalNumSeg
    temp <- f@exprs[(Segment*(k-1)+1):nrow(f@exprs), Time.loc]
    timeCentres[k] <- mean(temp)

    for ( j in CleanChan.loc){
        segSummary <- matrix(0, totalNumSeg, 8)
        # Calculating the eight features for each segment except the last one
        #  and store them in the matrix
        for ( k in 1:(totalNumSeg-1)){
            temp <- f@exprs[Segment*(k-1)+(1:Segment), c(j)]
            segSummary[k, ] <-
                c(
                quantile(temp, probs = c(5,20,50,80,95)/100),
                mean(temp),
                sapply(2:3, function(x){
                    moment(temp, order = x, center = TRUE)
                    })
                )
        }
        # Calculating the eight features for the last segment and store them
        #  in the matrix
        k <- totalNumSeg
        temp <- f@exprs[(Segment*(k-1)+1):nrow(f@exprs), c(j)]
        segSummary[k, ] <-
            c(
            quantile(temp, probs = c(5,20,50,80,95)/100),
            mean(temp),
            sapply(2:3, function(x){
                moment(temp, order = x, center = TRUE)
                })
            )

        storeMeans[[j]] <- segSummary[ ,6]
        quantiles[[j]] <- quantile(f@exprs[ ,j], probs = c(0,2,98,100)/100)
        meanRangePerc[j] <-
            (max(storeMeans[[j]]) - min(storeMeans[[j]])) /
            (quantiles[[j]]["98%"]-quantiles[[j]]["2%"])

        if (FirstOrSecond == "First"){
            segSummary <- rbind(sapply(1:ncol(segSummary), function(x) {
                (segSummary[ ,x]-mean(segSummary[ ,x])) / sd(segSummary[ ,x])
            } ))
            segSummary <- replace(segSummary, is.nan(segSummary), 0)
            cellDeleteExpo <- abs(segSummary)

            if (j==1){
                return(cellDeleteExpo)
            }

            # Sum up each row of the resulting matrix
            cellDeleteExpo <- sapply(1:nrow(cellDeleteExpo), function(x) {
                sum(cellDeleteExpo[x, ])
            } )
            cellDelete[[j]] <- cellDeleteExpo

            temp.vect <- rep(0, length(storeMeans[[j]])-1)
            for(l in 1:(length(storeMeans[[j]])-1) ){
                temp.vect[l] <-
                    abs(storeMeans[[j]][l]-storeMeans[[j]][[l+1]]) /
                    (quantiles[[j]]["98%"]-quantiles[[j]]["2%"])
            }
            maxDistJumped[j] <- max(temp.vect)
        }
    } # End of for-loop

    if (length(cellDelete) > 0){
        choosenChans <- NULL
        if (UseOnlyWorstChannels == TRUE){
            if (length(which(!is.na(meanRangePerc))) >= AmountMeanRangeKeep && AmountMeanRangeKeep >= 1){
                choosenChans <- sort(
                    unique(c(choosenChans, sapply(sort(meanRangePerc, decreasing = TRUE)[1:AmountMeanRangeKeep],
                        function(x) {which(x == meanRangePerc)} )))
                )
            }

            sdMeans <- sapply(1:length(storeMeans), function(x) {
                sd(storeMeans[[x]])
            } )

            if( length(which(!is.na(sdMeans))) >= AmountMeanSDKeep && AmountMeanSDKeep >= 1 ) {
                choosenChans <- sort(
                    unique(c(choosenChans, sapply(sort(sdMeans, decreasing = TRUE)[1:AmountMeanSDKeep],
                        function(x) {which(x == sdMeans)})))
                )
            }
            if ( Verbose == TRUE) {
                cat(paste0("The channels that are selected for cutting are: ",
                    paste0(f@parameters@data$desc[choosenChans],collapse=", "),
                        " (", paste0(choosenChans, collapse=","), ").\n")
                )
            }

            cellDelete <- cellDelete[choosenChans]
        }

        cellDelete1 <- do.call(cbind, cellDelete)

        cellDelete2 <- NULL
        for ( m in 1:length(cellDelete1[ ,1]) ){
            cellDelete2[m] <- sum(cellDelete1[m, ])
        }

        typeOfGating <- NULL
        # Check if the file is nice, before removing anything.
        # If it is, then we can avoid removing slivers.
        if ((round( max(maxDistJumped, na.rm=TRUE), digits=3) < MaxContin  )&&
            (round(mean(meanRangePerc, na.rm=TRUE), digits=3) < MeanOfMeans)&&
            (round( max(meanRangePerc, na.rm=TRUE), digits=3) < MaxOfMeans )&&
            is.null(GateLineForce)){
            densityGateLine <- NULL

            range_7SD <-  c(mean(cellDelete2)-7*sd(cellDelete2),
                            mean(cellDelete2)+7*sd(cellDelete2))
            densityGateLine <- range_7SD[2]
            deletedSegments <-    union(which(cellDelete2 < range_7SD[1]),
                                        which(cellDelete2 > range_7SD[2]))
            typeOfGating <- paste0(typeOfGating, "7SD")
            cellDelete2.org <- cellDelete2
        } else {
            # Finding outliers using by plotting density of the 8 measures,
            #  the index of cells that have 8 measures significantly different
            #  than the rest will be returned
            peaks_info  <- getPeaks(cellDelete2, tinypeak.removal = 0.1)
            peaks_xind  <- peaks_info$Peaks
            peaks_value <- which(density(cellDelete2)$x %in% peaks_xind)
            peaks_value <- density(cellDelete2)$y[peaks_value]
            peak_max    <- max(density(cellDelete2)$y)

            cellDelete2.org <- cellDelete2
            # Create very similar data that wont have a very large value in the
            #  distribution if there was one very obvious outlier.
            #  this allows deGate to be able to gate without resolution issues.
            cellDelete2 <- cellDelete2[which( cellDelete2 <=
                            (mean(cellDelete2)+5*sd(cellDelete2)) )]

            all_cut <- deGate(cellDelete2, tinypeak.removal = 0.001,
                            all.cuts = TRUE, upper = TRUE, verbose = F, count.lim = 3)

            if (!anyNA(all_cut)){
                # We need to set a limit on peaks
                #  because we don't want to overcut
                if (length(peaks_xind) > 1){
                    sn <- NULL
                    for ( ml in 1:(length(peaks_xind)-1) ){
                        ind <- intersect   (which(all_cut > peaks_xind[ml  ]),
                                            which(all_cut < peaks_xind[ml+1]))
                        y_ind <- which(density(cellDelete2)$x %in% all_cut[ind])
                        y_ind <- density(cellDelete2)$y[y_ind]

                        if (abs(peaks_value[ml+1]-min(y_ind)) <
                            abs(max(peaks_value)-min(y_ind)) * 0.015){
                            sn <- c(sn, ind)
                        }
                    }
                    # Remove the peaks that are too small to
                    #  be counted as a separate population
                    if (length(sn)>0){
                        peaks_xind <- peaks_xind[-(sn+1)]
                        all_cut <- all_cut[-sn]
                        # Remove peaks
                        typeOfGating <- paste0(typeOfGating, "RemPeaks")
                    }
                }
            }

            all_cut_value <- which(density(cellDelete2)$x %in% all_cut)
            all_cut_value <- density(cellDelete2)$y[all_cut_value]
            if (anyNA(peaks_xind)){
                # In case the density is a spike
                densityGateLine <- mean(density(cellDelete2)$x)/20
            } else {
                if (length(peaks_xind)==1){
                    upper_cut <- deGate(cellDelete2, tinypeak.removal = 0.1, upper = TRUE, use.upper = TRUE,
                                        alpha = 0.1, percentile = .95, verbose = F, count.lim = 3)
                    if (!is.na(upper_cut)){
                        # We don't want the upper cut and the all cut to be the
                        #  same because that means we are cutting at the 95 percentile
                        if (length(all_cut)==1 && upper_cut==all_cut){
                            upper_cut <- deGate(cellDelete2, tinypeak.removal = 0.1, upper = TRUE, use.upper = TRUE,
                                                alpha = 0.05, verbose = F, count.lim = 3)
                            if (upper_cut==all_cut){
                                upper_cut <- deGate(cellDelete2, tinypeak.removal = 0.1, upper = TRUE, use.upper = TRUE,
                                                    alpha = 0.01, verbose = F, count.lim = 3)
                            }
                        }
                    }

                    if (length(which(all_cut_value < MaxValleyHgt*peak_max))==0){
                        all_cut_min <- all_cut[1]
                    } else {
                        all_cut_min <- all_cut[min(which(all_cut_value < MaxValleyHgt*peak_max))]
                    }

                    if (!is.na(upper_cut) && !is.na(all_cut_min)){
                        if (upper_cut >= all_cut_min){
                            typeOfGating <- paste0(typeOfGating, "Upper")
                        } else {
                            typeOfGating <- paste0(typeOfGating, "AllCutMin")
                        }
                    }

                    densityGateLine <- max(upper_cut, all_cut_min, na.rm = TRUE)

                    if( densityGateLine == -Inf){ densityGateLine <- NA } # max returns -Inf if both were NA
                } else {
                    sind <- NULL
                    for ( t in all_cut ){ # Save the index that are less than the population cut off.
                        if (length(which(cellDelete2 > t)) <= MaxPercCut * length(cellDelete2)){
                            sind <- c(sind, t)
                        }
                    }
                    densityGateLine <- sind[1] # Use the min of the index for cut off
                    typeOfGating <- paste0(typeOfGating, "MaxPercCut")
                }
            }
            if (!is.null(GateLineForce)){
                densityGateLine <- GateLineForce
            }
            deletedSegments <- which(cellDelete2.org > densityGateLine)
        }
    }

    if (FirstOrSecond == "Second"){
        return(list(deletedSegments=0, quantiles=quantiles, storeMeans=storeMeans,
                    meanRangePerc=meanRangePerc, timeCentres=timeCentres))
    }
    deletedSegments <- suppressWarnings(sort(unique(deletedSegments)))

    if (length(deletedSegments) != 0) {
        # Split the consectutive sections into separate elements in the list
        range.seg.rem <- split(deletedSegments, cumsum(c(1, diff(deletedSegments) != 1)))

        size.in.a.row <- 5 # If greater or equal to 5 segments in a row, add 20% to each side
        adding.Segments <- NULL

        length.range.seg.rem <- sapply(1:length(range.seg.rem), function(x){ length(range.seg.rem[[x]]) })
        range.seg.rem.only.5.or.more <- range.seg.rem[which(length.range.seg.rem >= size.in.a.row)]

        range.seg.keep <- setdiff(1:totalNumSeg, unlist(range.seg.rem.only.5.or.more))
        range.seg.keep <- split(range.seg.keep, cumsum(c(1, diff(range.seg.keep) != 1)))


        if (length(range.seg.rem.only.5.or.more) >= 1){
            # Make sure range.seg.keep is first (will be skipped later because it is only length of 1)
            if (min(unlist(range.seg.rem.only.5.or.more)) <= min(unlist(range.seg.keep))){
                range.seg.keep <- c(0, range.seg.keep)
            }
            # Make sure range.seg.keep is last (will be skipped later because it is only length of 1)
            if (max(unlist(range.seg.rem.only.5.or.more)) >= max(unlist(range.seg.keep))){
                range.seg.keep <- c(range.seg.keep, max(unlist(range.seg.rem.only.5.or.more))+1)
            }
            # Adding buffer region only when the number of consecutive segments in the group to be deleted are quite large

            for ( b2 in 1:length(range.seg.rem.only.5.or.more) ){
                length_section <- length(range.seg.rem.only.5.or.more[[b2]])
                for (j1 in CleanChan.loc){

                    left.pop <- storeMeans[[j1]][range.seg.keep[[b2]]]
                    rght.pop <- storeMeans[[j1]][range.seg.keep[[b2+1]]]
                    rght.pop <- rev(rght.pop)

                    main.pop <- storeMeans[[j1]][range.seg.rem.only.5.or.more[[b2]]]

                    for ( side.pop in list(left.pop, rght.pop)){
                        if (length(side.pop) >= 2){
                            left.side <- identical(side.pop, left.pop)
                            # Code below works for good part to be lower than bad part, so we might need to transform our data
                            if ( mean(side.pop) >= mean(main.pop)){
                                side.pop <- -(side.pop-mean(main.pop))+mean(main.pop)
                            }

                            for (k1 in 1:length(side.pop)){
                                if (length(side.pop) <= 2 || sd(side.pop) == 0 || (max(side.pop)-min(side.pop) == 0) ){
                                    break;
                                }
                                if ( tail(side.pop, n=1) >= (mean(side.pop) + sd(side.pop)) ){
                                    if (left.side){
                                        adding.Segments <- c(adding.Segments, length(side.pop))
                                    } else {
                                        adding.Segments <- c(adding.Segments, range.seg.keep[[b2+1]][k1])
                                    }
                                    # Remove last point because we want to recalculate mean and sd on progressively cleaner data
                                    side.pop <- side.pop[-length(side.pop)]
                                } else {
                                    break;
                                }
                            }
                        }
                    }
                }
                adding.Segments <-
                    c(adding.Segments,
                        (min(range.seg.rem.only.5.or.more[[b2]])-ceiling(length_section/5)):
                          (min(range.seg.rem.only.5.or.more[[b2]])-1),
                        (max(range.seg.rem.only.5.or.more[[b2]])+1):
                          (max(range.seg.rem.only.5.or.more[[b2]])+ceiling(length_section/5))
                    )
            }
        }
        if (length(adding.Segments) >= 1){
            deletedSegments <- sort(unique(c(deletedSegments, adding.Segments[intersect(which(adding.Segments >=1),
                                                                                        which(adding.Segments <= totalNumSeg))])))
        }
    }
    return(list(deletedSegments=deletedSegments, quantiles=quantiles, storeMeans=storeMeans, typeOfGating=typeOfGating,
                meanRangePerc=meanRangePerc, timeCentres=timeCentres, densityGateLine=densityGateLine,
                cellDelete2=cellDelete2.org, choosenChans=choosenChans))
}


##################################################
# Time Function                                  #
##################################################
# Prints out the time since startTime
TimePrint <- function(startTime) {
    startTime <- as.POSIXct(startTime)
    difft <- difftime(Sys.time(), startTime, units="secs")
    format(.POSIXct(difft, tz="GMT"), "%H:%M:%S")
}