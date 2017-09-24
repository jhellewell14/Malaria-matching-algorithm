#' Reads in csv file data for humans
#'
#' @param fname string input filename
#' @param agefile string input file containing corresponding human ages
read_humans <- function(fname,agefile){

  # clear out rows without markers
  h <- read.table(fname,header=TRUE,sep=",") %>% filter(Marker!="" & Size>50)

  # add household ID and number of residents to samples
  h %<>% mutate(household=substr(Human.ID,1,3)) %>% group_by(household) %>% mutate(hnum=length(unique(Human.ID))) %>% ungroup


  # Match ages to humans
  # ages <- read.table("Copy of Age_AFIRM_2_Burkina.csv",header=TRUE,sep=",") %>% select(Human.ID=StudyID,Age)
  # h %<>% left_join(ages) %>% mutate(age_cat=ifelse(Age<5,"<5",ifelse(Age<15,"<15",">15")))

  return(h)
}



#' Reads in csv file of mosquito data
#'
#' @param fname string input of filename
read_mosquitoes <- function(fname){

  # Read in csv file and clean
  m <- read.table(fname,header=TRUE,sep=",",colClasses = c("character","character","character","numeric","numeric"))

  m %<>% filter(Marker!="" & Marker!="." & Height>50) %>% mutate(ID=paste(MosquitoID,Date.of.Collection,sep=""))

  m %<>% mutate(household=substr(as.character(ID),1,3))

  return(m)
}

#' Compares allele sizes at a single loci
#'
#' @param hm Human marker
#' @param mm Mosquito marker
compare_marker <- function(hm,mm){ # compares alleles at a single loci
  nh <- nrow(hm) # number of alleles in each sample
  nm <- nrow(mm)

  if(any(c(nh,nm)==0)){
    return(NA) # if no alleles at this loci, can't compare
  }

  if(any(c(nh,nm)==1)){
    dist <- (as.numeric(mm$Size) - hm$Size)^2
    return(min(dist))
  }

  te <- expand.grid(1:nm,1:nm) %>% filter(Var1!=Var2) # creates all the possible allele combinations
  te$dist1 <- (as.numeric(mm$Size[te$Var1]) - hm$Size[1])^2
  te$dist2 <- (as.numeric(mm$Size[te$Var2]) - hm$Size[2])^2
  te$dist <- (te$dist1 + te$dist2)/2
  return(min(te$dist))

}

#' Compares all 10 markers in a mosquito to all 10 markers in a human
#' @param human Table of marker sizes for human
#' @param mosquito Table of marker sizes for mosquito
compare_mosquito <- function(human,mosquito){ # compares a single mosquito and human
  m_names <- c("AMEL","D10S1248","D12S391","D19S433","D1S1656",
               "D22S1045","D2S1338","D2S441","D6S1043","TH01")

  dist_vec <- vapply(m_names,FUN=function(x){
    htab <- human %>% filter(Marker==x) # for each marker, selects allele sizes and
    mtab <- mosquito %>% filter(Marker==x) # sends to compare_marker
    return(compare_marker(htab,mtab))
  },FUN.VALUE = 1)
  av_dist <- mean(dist_vec,na.rm=TRUE) # finding the mean of the squared distances
  return(data.frame(dist=av_dist))
}

#' Compares all 17 markers in a mosquito to all 17 markers in a human
#'
#' @param human Table of human marker sizes
#' @param mosquito Table of mosquito marker sizes
compare_mosquito_16 <- function(human,mosquito){ # compares a single mosquito and human
  # m_names <- c("AMEL","D10S1248","D12S391","D19S433","D1S1656",
  #              "D22S1045","D2S1338","D2S441","D6S1043","TH01")
  m_names <- c("D10S1248","vWA","D16S539","D2S1338","AMEL","D8S1179","D21S11","D18S51",
               "D22S1045","D19S433","TH01","FGA","D2S441","D3S1358","D1S1656","D12S391","SE33")



  dist_vec <- vapply(m_names,FUN=function(x){
    htab <- human %>% filter(Marker==x) # for each marker, selects allele sizes and
    mtab <- mosquito %>% filter(Marker==x) # sends to compare_marker
    return(compare_marker(htab,mtab))
  },FUN.VALUE = 1)
  av_dist <- mean(dist_vec,na.rm=TRUE) # finding the mean of the squared distances
  return(data.frame(dist=av_dist))
}

#' Compares a mosquito to all humans that live in the house within which it was caught
#'
#' @param mosquito Mosquito marker data
#' @param humans All human marker data
#' @param threshold Distance below which a mosquito and human are declared a match
#' @param sixteen boolean value denoting whether 10 or 17 markers are being used
human_matches <- function(mosquito,humans,threshold,sixteen){ # matches a mosquito to humans from same household
  hum_house <- unique(mosquito$household)
  rel_hum <- humans %>% dplyr::filter(household==hum_house) # select human samples in same house

  # If no humans found in the same household, search all humans
  if(nrow(rel_hum)==0){rel_hum <- humans}

  # are there >2 alleles at >2 loci?
  comp_flag <- mosquito %>% group_by(Marker) %>% count %>% mutate(comp=ifelse(n>2,TRUE,FALSE))  %>% .$comp %>% sum()>2

  # are there enough alleles to make a comparison?
  if(mosquito %>% .$Marker %>% unique %>% length<6){return(data.frame("verdict"="No Amplification",distance=NA,multiple=FALSE))}

  # find allele combinations with min distance for each human
  if(sixteen==FALSE){
    rel_hum %<>% group_by(Human.ID) %>% do(compare_mosquito(human=.,mosquito = mosquito))
  }else{
    rel_hum %<>% group_by(Human.ID) %>% do(compare_mosquito_16(human=.,mosquito = mosquito))
  }

  # depending on if mosquito sample seems complex or not, return 1 match or several
  if(comp_flag==FALSE){ # simple case, return best match
    min_dist <- rel_hum %>% ungroup %>% top_n(-1,wt=dist)
    if(min_dist$dist <= threshold){
      verdict <- min_dist$Human.ID
      return(data.frame("verdict"=verdict,"distance"=min_dist$dist,"multiple"=comp_flag))
    }else{
      verdict <- "Outside"
      return(data.frame("verdict"=verdict,"distance"=min_dist$dist,"multiple"=comp_flag))
    }
  }else{ # complex case, return all matches under threshold
    # if(any(rel_hum$dist<=threshold)==FALSE){ # no matches found under threshold
    #   min_dist <- rel_hum %>% ungroup %>% top_n(-1,wt=dist)
    #   #verdict <- "Outside"
    #   verdict <- min_dist$Human.ID
    #   return(data.frame("verdict"=verdict,"distance"=min_dist$dist,"multiple"=comp_flag))
    # }else{ # matches found under threshold
    #   rel_hum %<>% ungroup %>% filter(dist<=threshold)
    #   return(data.frame("verdict"=rel_hum$Human.ID,"distance"=rel_hum$dist,"multiple"=rep(comp_flag,nrow(rel_hum))))
    # }
    return(data.frame("verdict"=rel_hum$Human.ID,"distance"=rel_hum$dist,"multiple"=rep(TRUE,nrow(rel_hum)))) #can be used to display
    # distances between a multiple feed mosquito and ALL humans in household (Bronner requested for sensitivity analysis)

  }
}


# human_matches_2 <- function(mosquito,humans,threshold){
#   hum_house <- unique(mosquito$household)
#   rel_hum <- humans %>% dplyr::filter(household==hum_house) # select human samples in same house
#   if(nrow(rel_hum)==0){rel_hum <- humans}
#
#   # select the two allele reads with the largest height at each loci
#   mosquito %<>% group_by(MosquitoID,Marker) %>% top_n(2,wt=Height) %>% ungroup()
#
#   # find minimum distances for each human
#   rel_hum %<>% group_by(Human.ID) %>% do(compare_mosquito(human=.,mosquito = mosquito))
#
#   # go through process for simple cases
#   min_dist <- rel_hum %>% ungroup %>% top_n(-1,wt=dist)
#   if(min_dist$dist <= threshold){
#     verdict <- min_dist$Human.ID
#     return(data.frame("verdict2"=verdict,"distance2"=min_dist$dist))
#   }else{
#     verdict <- "Outside"
#     return(data.frame("verdict2"=verdict,"distance2"=min_dist$dist))
#   }
# }

#' Performs whole analysis for given mosquito and human datasets
#'
#' @param mos All mosquito marker sizes
#' @param hum All human marker sizes
#' @param threshold Maximum allowable distance for a match
#' @param sixteen Boolean denoting 10 or 17 markers used
compare_all <- function(mos,hum,threshold,sixteen){

  # First type of matching, return all matches found for suspected multiple feeds
  res <- mos %>% group_by(ID) %>% do(human_matches(mosquito=.,humans=h,threshold=threshold,sixteen))

  # Second type of matching, for multiple feeds remove the allele at each loci with the lowest heights (potential artifacts)
  # and match as a simple case

  # re_run <- unique(res$ID[which(res$multiple==TRUE)])
  # res2 <- mos %>% filter(ID %in% re_run) %>% group_by(ID) %>% do(human_matches_2(mosquito=.,humans=h,threshold=threshold))
  #
  # combine results
  #out <- left_join(res,res2)
  out <- res

  return(out)
}

# unique_at_locus <- function(temp){
#     sz <- temp$Size
#     k <- nrow(temp)
#     n <- temp$hnum[1]*2
#     temp %<>% mutate(ver=jd(Size,sizes=sz)) %>% filter(ver==TRUE)
#     if(k >= n){
#       return(temp)
#     }else{
#       return(temp %>% filter(age_cat=="No lol"))
#     }
# }



  # x <- h %>% filter(household=="29C")
  # marktab <- x %>% filter(Marker=="D10S1248" & !is.na(age_cat))
  # marktab$verdict <- rep(FALSE,nrow(marktab))
  # lw <- marktab %>% filter(age_cat=="<5") %>% .$Size
  # bt <- marktab %>% filter(age_cat=="<15") %>% .$Size
  # hg <- marktab %>% filter(age_cat==">15") %>% .$Size
  #
  # for(i in 1:nrow(marktab)){
  #   if(marktab$age_cat[i]==">5"){
  #     if(all(abs(marktab$Size[i] - bt)>2) & all(abs(marktab$Size[i] - hg)>2)) marktab$verdict[i] >- TRUE
  #   }
  #   if(marktab$age_cat[i]==">15"){
  #     if(all(abs(marktab$Size[i] - lw)>2) & all(abs(marktab$Size[i] - hg)>2)) marktab$verdict[i] >- TRUE
  #   }
  #   if(marktab$age_cat[i]==">15"){
  #     if(all(abs(marktab$Size[i] - bt)>2) & all(abs(marktab$Size[i] - lw)>2)) marktab$verdict[i] >- TRUE
  #   }
  # }
  # marktab



# unique_age_cat_household <- function(tab){
#   tab %>% group_by(age_cat)
# }


# jd <- function(val,sizes){
#   ret <- vapply(val,FUN.VALUE = TRUE,FUN = function(x){
#     ifelse(sum(abs(x-sizes)<2)==1,TRUE,FALSE)
#   })
#   return(ret)
# }

# check_household_mosquitoes <- function(human,mosquitoes){
#   mosquitoes %<>% filter(household==human$household & Marker == human$Marker & Size == human$Size)
#   if(nrow(mosquitoes)>1){
#     ret <- data_frame(Human.ID=human$Human.ID,Marker=human$Marker,age_cat=human$age_cat,Mosquito.ID=mosquitoes$MosquitoID)
#   }else{
#     ret <- data.frame()
#   }
#   return(ret)
# }