#load packages
library('tidyverse')


### FUNCTION TO CREATE AN EMPTY PSD PLOT
f_PSD_empty <- function(
  height_labelbar = 0.075, #height of single label bar above plot
  percentage = T,          #plot percentrages rather than fractions
  box_thickness = 0.25,    #thickness of label boxes
  textsize_type = 2.5,     #size of labels, major type
  textsize_subtype = 2,    #size of labels, minor type
  xlim = NA
){

  #Particle size distribution classes
  dc <- tibble(
    type = c('CLAY','SILT','SILT','SILT','SAND','SAND','SAND','GRAVEL','GRAVEL','GRAVEL','COBBLES','BOULDERS'),
    subtype = c(NA,'Fine','Medium','Coarse','Fine','Medium','Coarse','Fine','Medium','Coarse',NA,NA),
    dmin = c(4e-05,0.002,0.006,0.02,0.06,0.2,0.6,2,6,20,60,200),
    dmax = c(0.002,0.006,0.02,0.06,0.2,0.6,2,6,20,60,200,2e+03)
  )
  #limit based on input
  if (any(is.na(xlim))){
    xlim <- c(min(dc$dmin),max(dc$dmax))
  } else {
    dc <- dc[dc$dmax>xlim[1] & dc$dmin<xlim[2],]
    dc$dmin <- pmax(dc$dmin, xlim[1])
    dc$dmax <- pmin(dc$dmax, xlim[2])
  }

  #function to get all minor breaks in log10_range
  f_log10_minorbreaks <- function(x){
    minx <- floor(min(log10(x), na.rm=T))-1
    maxx <- ceiling(max(log10(x), na.rm=T))+1
    n_major <- maxx-minx+1
    major_breaks <- seq(minx,maxx,by=1)
    minor_breaks <- 
      rep(log10(seq(1,9,by=1)), times=n_major)+
      rep(major_breaks,each=9)
    return(10^(minor_breaks))
  }

  #add plot data to classes - major classes labels
  dc1 <- dc %>%
    group_by(type) %>%
    summarise(
      dmin = min(dmin),
      dmax = max(dmax),
      ymin = ifelse(any(is.na(subtype)),1,1+height_labelbar)
    ) %>% 
    mutate(
      ymax = 1+2*height_labelbar,
      dlabel = 10^(0.5*(log10(dmin)+log10(dmax))),
      ylabel = 0.5*(ymin+ymax)
    )
  #add plot data to classes - minor classes labels
  dc2 <- dc[!is.na(dc$subtype),] %>%
    group_by(type,subtype) %>%
    summarise(
      dmin = min(dmin),
      dmax = max(dmax)
    ) %>% 
    mutate(
      ymin = 1,
      ymax = 1+1*height_labelbar,
      dlabel = 10^(0.5*(log10(dmin)+log10(dmax))),
      ylabel = 0.5*(ymin+ymax)
    )

  #Create multiplier if percentages are used
  if (percentage==T){
    mult <- 100
  } else {
    mult <- 1
  }
  
  #plot empty plot
  p <- ggplot() + 
    theme_bw() +
    geom_rect(data=dc1, aes(xmin=dmin, xmax=dmax, ymin=mult*ymin, ymax=mult*ymax), fill='white', color='black', size=box_thickness) +
    geom_rect(data=dc2, aes(xmin=dmin, xmax=dmax, ymin=mult*ymin, ymax=mult*ymax), fill='white', color='black', size=box_thickness) +
    geom_text(data=dc1, aes(x=dlabel, y=mult*ylabel, label=type), hjust=0.5, vjust=0.5, size=textsize_type) +
    geom_text(data=dc2, aes(x=dlabel, y=mult*ylabel, label=subtype), hjust=0.5, vjust=0.5, size=textsize_subtype) +
    scale_y_continuous(
      limits=c(0,mult*(1+2*height_labelbar)), 
      breaks=mult*seq(0,1,0.2),
      expand=c(0,0)
    ) + 
    scale_x_log10(
      limits=xlim,
      breaks=scales::trans_breaks("log10", function(x) 10^x), 
      labels=scales::trans_format("log10", scales::math_format(10^.x)),
      minor_breaks=f_log10_minorbreaks(xlim),
      expand=c(0,0)
    ) + 
    annotation_logticks(side="b")
  
  #add labels
  p <- p + xlab('Particle size [mm]')
  if (percentage==T){
    p <- p + ylab('Percentage passing [%]')
  } else {
    p <- p + ylab('Fraction passing [-]')
  }
  
  #return
  return(p)
}


### FUNCTION TO ADD DATA TO PLOT
f_PSD_add_data <- function(p, d,
  legend=F,
  percentage=T,
  markers=T,
  palette='Set1')
{
  #Create multiplier if percentages are used
  if (percentage==T){
    mult <- 100
    d$y <- mult*d$y
  } else {
    mult <- 1
  }
  #traces: sort by first field that is not <x> or <d>
  col_sort <- setdiff(colnames(d),c('d','y'))
  #add traces
  if (length(col_sort)==0){
    if (markers==T){
      p <- p + geom_point(data=d, aes(x=d, y=y, color='A'))
    }
    p <- p + geom_line(data=d, aes(x=d, y=y, color='A', linetype='A'))

  } else {
    if (markers==T){
      p <- p + geom_point(data=d, aes_string(x='d', y='y', color=col_sort[1]))
    }
    p <- p + 
      geom_line(data=d, aes_string(x='d', y='y', color=col_sort[1], linetype=col_sort[1]))
  }
  p <- p + 
    scale_color_brewer(name=legend, palette=palette) +
    scale_linetype_discrete(name=legend)
  #suppress legend
  if (legend==F){
    p <- p + theme(legend.position='none')
  }
  #return
  return(p)
}

#Add D10 etc arrows
f_PSD_add_D10 <- function(p,  #empty psd plot
                          d,  #dataframe with fields <d> for diameter and <y> for psd (in fraction)
                          yi, #array with interpolation values (0.10, 0.30 etc)
                          nsignif=2, 
                          percentage=T, 
                          arrowlength=0.2,
                          label_size=4,
                          xlim=NA){
  #limit based on input
  if (any(is.na(xlim))){
    xlim <- c(4e-05,2e3)
  }
  #Create multiplier if percentages are used
  if (percentage==T){
    mult <- 100
  } else {
    mult <- 1
  }
  #interpolate data to find diameters
  di <- data.frame(
    y = mult*yi,
    d = 10^approx(x=d$y, y=log10(d$d), xout=yi, method='linear')$y
  )
  di$label <- paste('~D[',yi*100, ']==', signif(di$d,nsignif), '~"mm"', sep='')
  #add some arrows to a plot
  for (i in 1:nrow(di)){
    p <- p + 
      annotate('segment', x=xlim[1],y=di$y[i],xend=di$d[i],yend=di$y[i], arrow=arrow(length=unit(arrowlength,"cm"))) +
      annotate('segment', x=di$d[i],y=di$y[i],xend=di$d[i],yend=0, arrow=arrow(length=unit(arrowlength,"cm"))) +
      annotate('text', x=di$d[i], y=di$y[i], hjust=0, vjust=0.5, label=di$label[i], parse=T, size=label_size)
    }
  #return
  return(p)
}

# FUNCTION TO INDICATE FRACTIONS IN PSD
f_PSD_indicatefractions <- function(
  p,  #current PSD plot
  d,   #dataframe with PSD data. fields diameter <d> and cumulative fraction <y>
  percentage=T,
  line_color='blue',
  line_type=2,
  line_size=0.5
){
  #percentage settings
  if (percentage==T){
    mult <- 100
  } else {
    mult <- 1
  }
  #lims
  lims <- c(0.002,0.063,2,63,200)
  names <- c('clay','silt','sand','gravel','cobbles','boulders')
  fracs <- rep(0,length(names))
  #interpolate
  y <- approx(log10(d$d),d$y,log10(lims), rule=1)$y
  fracs <- diff(c(0,y,1))
  df <- setNames(data.frame(matrix(fracs,nrow=1)),names)
  #add lines in PSD
  p <- p + 
    geom_vline(xintercept=lims, color=line_color, linetype=line_type, size=line_size) +
    geom_hline(yintercept=mult*y, color=line_color, linetype=line_type, size=line_size) 
  #return
  return(p)
}

# FUNCTION TO CALCULATE MASS FRACTIONS FROM PSD
f_PSD_getmassfractions <- function(
  d  #dataframe with PSD data. fields diameter <d> and cumulative fraction <y>
)
{
  #traces: sort by first field that is not <x> or <d>
  col_sort <- setdiff(colnames(d),c('d','y'))
  #particle size limits
  lims <- c(0.002,0.063,2,63,200,630)
  names <- c('clay','silt','sand','gravel','cobble','boulder','large boulder')
  fracs <- rep(0,length(names))
  #interpolate each trace at specified diameters
  d %>%
    group_by(.dots=col_sort) %>%
    group_modify(~data.frame(d=lims,y=approx(log10(.x$d), .x$y, lims)$y))
  
  
  y <- approx(log10(d$d),d$y,log10(lims), rule=1)$y
  fracs <- diff(c(0,y,1))
  #incomplete data 
  #create dataframe
  return(setNames(data.frame(matrix(fracs,nrow=1)),names))
}


#test
#d <- tibble(
#  d=c(10^seq(-3,2,l=51), 10^seq(-3,2,l=51)),
#  y=c(0.5+0.5*tanh(log10(10^seq(-3,2,l=51))/1), 0.5+0.5*tanh(log10(10^seq(-3,2,l=51))/10)),
#  l=c(rep('A',51),rep('B',51))
#)
#
#p <- f_PSD_empty()
#p1 <- f_PSD_add_data(p, d, marker=F)
#p2 <- f_PSD_indicatefractions(p1, d)
#p2
