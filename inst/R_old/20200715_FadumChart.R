library('tidyverse')
library("RColorBrewer")

f_fadumschart_empty <- function(
  val_inf = 1e6, #value for infinite approximation
  palette='Dark2'
){

  #colors
  #colo <- brewer.pal(n=4, name=palette) #colors
  colo <- c('black','darkred','darkblue','darkgreen')
  
  #specify values
  d0 <- data.frame(yz=c(0.5,1.0,val_inf), x0=0.01, x1=10, color=colo[1])
  d1 <- data.frame(yz=c(0.1,0.2,0.3,0.4), x0=0.03, x1=10, color=colo[2])
  d2 <- data.frame(yz=c(0.6,0.7,0.8,0.9,1.2), x0=0.1, x1=10, color=colo[3])
  d3 <- data.frame(yz=c(1.4,1.6,1.8,2.0,2.5,3.0), x0=1, x1=10, color=colo[4])
  d <- rbind(d0,d1,d2,d3)
  d <- d[order(d$yz),]
  
  #expand - long data
  dl <- rbind(
    expand.grid(yz=d0$yz, xz=10^seq(log10(d0$x0[1]),log10(d0$x1[1]),l=101)),
    expand.grid(yz=d1$yz, xz=10^seq(log10(d1$x0[1]),log10(d1$x1[1]),l=101)),
    expand.grid(yz=d2$yz, xz=10^seq(log10(d2$x0[1]),log10(d2$x1[1]),l=101)),
    expand.grid(yz=d3$yz, xz=10^seq(log10(d3$x0[1]),log10(d3$x1[1]),l=101))
  )
  
  #calculation functions
  fI <- function(xz,yz){
    return(1/(2*pi)*(xz*yz*(xz^2+yz^2+2)/((1+xz^2)*(1+yz^2)*sqrt(1+xz^2+yz^2))+atan(xz*yz/sqrt(1+xz^2+yz^2))))
  }
  
  #cals
  dl$I <- fI(dl$xz, dl$yz)
  
  #labels
  xmin <- min(d$x0)
  xmax <- max(d$x1)
  marg <- 0.06
  d$xlab <- 10^(log10(xmin) + (1.0-marg)*(log10(xmax)-log10(xmin)))
  d$ylab <- fI(d$xlab, d$yz)
  d$text <- d$yz
  d$text[d$yz==val_inf] <- 'infinity'
  d$text[d$yz==min(d$yz)] <- paste('B/z==',min(d$yz),sep='')
  #offset
  offset <- 0.04
  noffset <- 7
  for (i in seq(noffset)){
    d$xlab[length(d$yz)-i+1] <- 10^(log10(xmin) + (1.0-marg-offset*(noffset-i))*(log10(xmax)-log10(xmin)))
  }
  d$ylab <- fI(d$xlab, d$yz)
  
  #plot
  p <- ggplot(dl, aes(x=xz, y=I, color=as.factor(yz))) + 
    theme_bw() + theme(legend.position='none') +
    geom_line(size=0.3, show.legend=F) +
    xlab("L/z [-]") +
    ylab("I [-]") +
    coord_cartesian(xlim=c(min(dl$xz),max(dl$xz)), ylim=c(0,0.26), expand=F) +
    scale_y_continuous(breaks=seq(0,0.26,by=0.02), 
                       minor_breaks=seq(0,0.26,by=0.01)) +
    scale_x_log10(breaks=c(0.01,0.1,1,10),
                  minor_breaks=c(0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,2,3,4,5,6,7,8,9)) +
    scale_color_manual(values=as.character(d$color)) + 
    geom_label(data=d,aes(x=xlab,y=ylab,label=text), hjust=0.5, vjust=0.5, size=2, label.padding=unit(0.1,"lines"), parse=T, show.legend=F)

  #return
  return(p)
}

f_fadumschart_adddata <- function(
  p, # fadums chart plotted
  dd,  # dataframe, with <Lz>, <Bz> and <label> (optional)
  legpos = NA,
  crosshairs = T
){
  #add labels
  if (!'label'%in%colnames(dd)){
    dd$label <- as.factor(seq(nrow(dd)))
  }
  #calculate I
  fI <- function(xz,yz){
    return(1/(2*pi)*(xz*yz*(xz^2+yz^2+2)/((1+xz^2)*(1+yz^2)*sqrt(1+xz^2+yz^2))+atan(xz*yz/sqrt(1+xz^2+yz^2))))
  }
  dd$I <- fI(dd$Lz, dd$Bz)
  #add crosshairs
  if (crosshairs==T){
    ddl <- dd %>%
      group_by(label) %>%
      group_modify(~data.frame(Lz=c(.x$Lz,.x$Lz,1e-6),I=c(0,.x$I,.x$I)))
    p2 <- p + 
      geom_path(data=ddl, aes(x=Lz,y=I,linetype=as.factor(label)), color='black')
  }
  #add data points 
  p2 <- p2 + 
    geom_point(data=dd, aes(x=Lz, y=I, shape=as.factor(label)), color='black')
  #legend
  if (!is.na(legpos)){
    p2 <- p2 + 
      theme(legend.position=legpos) + 
      scale_linetype_discrete(name='') + 
      scale_shape_discrete(name='')
  }
  #return
  return(p2)
}

#Test
 #d <- data.frame(Lz=c(1,2),Bz=c(1,2))
 #p <- f_fadumschart_empty()
 #p2 <- f_fadumschart_adddata(p,d,legpos='right')
 #p2
