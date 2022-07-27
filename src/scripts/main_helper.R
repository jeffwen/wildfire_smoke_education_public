library(purrr)

######### general helpers ######### 
## from https://stackoverflow.com/a/56553052
mode <- function(x) {
  ux <- unique(na.omit(x))
  tx <- tabulate(match(x, ux))
  if(length(ux) != 1 & sum(max(tx) == tx) > 1) {
    if (is.character(ux)) return(NA_character_) else return(NA_real_)
  }
  max_tx <- tx == max(tx)
  return(ux[max_tx])
}

## similar to %in%
`%not_in%` <- purrr::negate(`%in%`)

######### Figure 1 ######### 
## helper for generating labels for figure 1 maps
dataPercentileLabeler <- function(levels_vector, round_digits=0) {
  
  new_levels <- imap_chr(levels_vector, function(., idx) {
    
    if (idx==1) {
      
      return(paste0("\u2264 ", round(as.double(str_split(str_sub(levels_vector[idx], 2, -2), ",")[[1]]), round_digits)[-1]))
      
    } else if (idx==length(levels_vector)){ 
      ## for the max element
      int_vect <- round(as.double(str_split(str_sub(levels_vector[idx], 2, -2), ",")[[1]]), round_digits)
      temp_val <- int_vect[1]+(1/10)^getDecimals(int_vect[1])
      
      ## need to figure out this...
      # if both are 0 then just return 0
      if (round(temp_val, 10)==0 | all(int_vect==0)){
        ## iterating to the next value so buckets dont overlap
        return(paste0(0))
      } else if(int_vect[-1]==0) {
        
        return(paste0(format(temp_val, scientific=FALSE)," - ", int_vect[-1])) 
      } else if(int_vect[1]==0){
        return(paste0(format((1/10)^round_digits, scientific=FALSE)," +"))
      } else{
        ## iterating to the next value so buckets dont overlap
        return(paste0(format(temp_val, scientific=FALSE)," +")) 
      }
      
    } else {
      ## for other entries based on the cut add one to the first element 
      int_vect <- round(as.double(str_split(str_sub(levels_vector[idx], 2, -2), ",")[[1]]), round_digits)
      
      lower_bound <- format(int_vect[1]+(1/10)^getDecimals(int_vect[1]), scientific=FALSE)
      upper_bound <- format(int_vect[-1], scientific=FALSE)
      
      return(paste(c(lower_bound,upper_bound), collapse=" - "))
    }
  })
  
  return(new_levels)
}

## get decimal place helper for map legend
## from https://r.789695.n4.nabble.com/how-to-calculate-a-numeric-s-digits-count-td4698742.html
getDecimals<-function(x)
{
  str<-format(x, scientific=FALSE)
  if(is.na(strsplit(str,"\\.")[[1]][2])) return(0)
  else return(nchar(strsplit(str,"\\.")[[1]][2]))  
}

######### Figure 2a ######### 
## bootstrap for avg with parallel lapply and fixest
getBootstrapAvg <- function(data, n=100, x="smoke_days", ctrls="avg_temp + avg_temp_2", 
                            fes="GEOID + year_grade", cluster="county_fips"){
  
  ## run the simulations
  coef_sim_list <- pbmcapply::pbmclapply(1:n, function(i){
    
    ## sample data by cluster
    newdata <- blockResample(data=data, cluster=cluster)
    
    ## create models as a list (just the coefficients)
    mod_coef_list <- modelGeneratorFixestAvg(data=newdata, x=x, ctrls=ctrls, fes=fes, cluster=cluster) 
    
    return(mod_coef_list)
  })
  
  ## combine the list of results by subject
  ## providing named list to help with named output list
  coef_list <- lapply(list("avg"="avg"), function(subject) {
    
    ## bind separate list into matrix
    do.call(rbind, lapply(coef_sim_list, function(coef_vec) 
    {coef_vec[[subject]]}))
  })
  
  return(coef_list)
}

## block resampling for bootstrap standard error/ ci estimation
blockResample <- function(data, cluster){
  
  ## sample by cluster id
  cluster_ids <- unique(data[[cluster]])
  block_idx <- base::sample(cluster_ids, size=length(cluster_ids), replace=T)
  
  ## merge back the observations that belong to each cluster
  block_dt <- data.table::data.table(col_name=block_idx)
  names(block_dt) <- cluster ## rename...probably a better way
  resampled_df <- data.table::merge.data.table(block_dt, data, by=cluster, allow.cartesian=T) %>%
    as.data.frame()
  
  return(resampled_df)
}

## model generator for avg of subjects
modelGeneratorFixestAvg <- function(data, x, ctrls, fes, cluster, weight) {
  
  ## define the formula for the model
  rhs_fmla <- paste0(x,'+', ctrls, "|", fes)
  
  ## create models
  avg_mod_res <- fixest::feols(as.formula(paste0("avg_score_unweighted",'~',rhs_fmla)), cluster=cluster, weights=data$avg_student_population,data=data)
  
  return(list("avg"=coef(avg_mod_res)))
}

## generate estimates for bootstrap over the support 
genRangeEstimateNew <-  function(clist, support_range, newdata, mid_idx=MID_IDX){
  
  ## for each subject
  result_dt <- lapply(seq_along(clist), function(i){
    
    ## estimate outcome for newdata values
    results <- clist[[i]] %*% t(newdata)
    
    ## center the data to midpoint
    results <- results - results[,mid_idx]
    
    ## convert to datatable with support
    results <- as.data.table(t(results))
    results$x <- support_range
    
    ## change from wide to long dataset
    results <- data.table::melt(results, id.vars=c("x"), variable.name="run", value.name="y")
    results$subject <- names(clist)[[i]]
    
    return(results)
  }) %>% rbindlist()
  
  return(result_dt)
}


######### Figure 2c ######### 
## random noise permutation test, shuffle obs within county randomly
permTest <- function(i, outcome, rhs, data, treatment, cluster_err, sample_group=NULL){
  
  ## sample and create new data with shuffled treatment
  
  if(!is.null(sample_group)){
    new_data <- data %>% 
      group_by(!!sample_group) %>% 
      mutate(!!paste0(treatment,"_new") := sample(get(c(treatment)), size=n())) %>% 
      ungroup()
    
  } else{
    new_data <- data %>% 
      mutate(!!paste0(treatment,"_new") := sample(get(c(treatment)), size=n())) %>% 
      ungroup()
  }
  
  fixest_df <- lapply(outcome, function(x){
    fmla <- paste0(x,"~",rhs)
    fixe_est <- feols(as.formula(stringr::str_replace_all(fmla, treatment, paste0(treatment,"_new"))),
                      weights=new_data$avg_student_population,
                      data=new_data)
    
    ## no clustering SEs here because we are just interested in the coef estimates
    data.frame("coef"=c(summary(fixe_est)$coefficients[[paste0(treatment,"_new")]]),
               "subject"=c(as.character(fixe_est$fml)[2])) %>% 
      mutate(subject=ifelse(subject=="avg_score_unweighted",
                            "avg",
                            ifelse(x=="mn_all_ela",
                                   "ELA",
                                   "Math")))
  }) %>%
    bind_rows()
  
  return(fixest_df)
}

######### Figure 4a ######### 
## get the coefficient of interest for categorical/categorical/continuous interaction
getcoef <- function(coef, mod1, mod2, treat_name, mod1_level, mod2_level,variable='Estimate'){
  return(coef[paste0(mod1,mod1_level,":",mod2,mod2_level,":",treat_name),variable])
}

## map plotting helper function
mapPlot <- function(plot_df, plot_col, fname, legend_title, color_palette=NA, outline_col=NA, 
                    width=1976, height=1142, num_classes=6, font_scaling=0.7, default_color="#737373", na_color="gray90", round_digits=4){
  
  ## determine if color palette was supplied
  if (any(is.na(color_palette))){
    ## viridis
    COLOR_PALETTE <- viridisLite::magma(n=num_classes, alpha=0.9, begin=0.1, end=0.9, direction=-1) 
  } else{
    ## other color palette
    COLOR_PALETTE <- color_palette
  }
  
  par(mar = c(2, 2, 2, 2), # Dist' from plot to side of page
      mgp = c(1, 0.4, 0), # Dist' plot to label
      ps = 14, cex = 1*FONT_SCALING, cex.main = 1*FONT_SCALING, # set the default font size to 12
      las = 1, # Rotate y-axis text
      tck = -.01, # Reduce tick length
      xaxs = "i", yaxs = "i",
      col=DEFAULT_COLOR) # Remove plot padding
  
  png(file=str_glue("{str_split(fname, '.pdf')[[1]][1]}_temp.png"), width=width, height=height)
  
  ## if outlines should be plotted
  if (!any(is.na(outline_col))){
    
    ## plot the NAs
    plot(plot_df[is.na(plot_df[plot_col]),][plot_col],
         col=na_color,
         border=na_color, lwd=0.02, main="", reset=FALSE)
    
    ## main plot
    plot(plot_df[plot_col],
         col=(COLOR_PALETTE)[plot_df[[plot_col]]],
         border="gray90", lwd=0.02, main="", add=TRUE)
    
    dev.off()
    
    
    ## create PDF from bitmap and legend as vector
    grDevices::cairo_pdf(file=fname,width=9.5,height=5.5)
    
    plot.new()
    
    grid.draw(rasterGrob(readPNG(str_glue("{str_split(fname, '.pdf')[[1]][1]}_temp.png"), native = FALSE),interpolate = FALSE))
    
    ## district/smoke legend
    l=legend("bottomleft",
             legend = c("NA", dataPercentileLabeler(levels(plot_df[[plot_col]]), round_digits=round_digits)),
             text.col = default_color,
             fill = c(na_color, COLOR_PALETTE),
             title=legend_title,
             title.adj=c(0),
             title.col=default_color,
             border=NA,
             bty = "n", # turn off the legend border
             cex = 1.1*font_scaling,
             inset=c(0.03,0),
             xjust=1,
             yjust=0
    ) # decrease the font / legend size
    
    ## outline
    plot(outline_col, col=alpha("black", 0), border="black", lwd=1, main="", add=TRUE)
    
  } else {
    ## plot the NAs
    plot(plot_df[is.na(plot_df[plot_col]),][plot_col],
         col=na_color,
         border=na_color, lwd=0.02, main="", reset=FALSE)
    
    ## main plot
    plot(plot_df[plot_col],
         col=(COLOR_PALETTE)[plot_df[[plot_col]]],
         border="gray90", lwd=0.02, main="", add=TRUE)
    
    dev.off()
    
    ## create PDF from bitmap and legend as vector
    grDevices::cairo_pdf(file=fname,width=9.5,height=5.5, pointsize=14)
    
    plot.new()
    
    grid.draw(rasterGrob(readPNG(str_glue("{str_split(fname, '.pdf')[[1]][1]}_temp.png"), native = FALSE),interpolate = FALSE))
    
    if(legend_title!=""){
      ## district/smoke legend
      l=legend("bottomleft",
               legend = c("NA", dataPercentileLabeler(levels(plot_df[[plot_col]]), round_digits=2)),
               text.col = DEFAULT_COLOR,
               fill = c(NA_COLOR, COLOR_PALETTE),
               border=NA,
               bty = "n", # turn off the legend border
               cex = 1.1*FONT_SCALING,
               inset=c(-0.07,-0.35),
               xpd=TRUE,
               xjust=1,
               yjust=0
      ) # decrease the font / legend size
      text(x=l$rect$left*(1/1.2), y=l$rect$top*1.48, 
           "Average test\nscore (stdev)",
           col=DEFAULT_COLOR,
           xpd=TRUE,
           adj=0,
           cex=1.1*FONT_SCALING)
       
    }
  }
  
  dev.off()
}


