library('tidyverse')
library("RColorBrewer")
library('directlabels')

#empty plot
f_trianglestresschart_empty <- function(
  xlim = c(-1,2),
  ylim = c(0,4),
  gridsize = 0.05,
  showgrid = T,
  val_inf = 1e6,   #value for infinite approximation,
  larrow = 0.3, #plot arrow length
  lsize = 3.5,
  labelsize = 3
){
  
  #colors
  colo <- brewer.pal(n=4, name='Dark2')
  colo <- c('black','darkred','darkblue','darkgreen')
  
  #specify values
  xB <- seq(xlim[1]-0.5,xlim[2]+0.5, gridsize)
  zB <- seq(ylim[1], ylim[2]+0.5, gridsize)
  
  #expand - long data
  dl <- expand.grid(xB=xB, zB=zB, KEEP.OUT.ATTRS=F)
  
  #calculation functions
  fI <- function(xB,zB){
    delta <- atan2(xB-1,zB)
    alphadelta <- atan2(xB,zB)
    alpha <- alphadelta - delta
    I <- (xB*alpha - 0.5*sin(2*delta))/pi
    return(I)
  }
  
  #cals
  dl$I <- fI(dl$xB,dl$zB)
  
  #labels
  xint <- 0.9
  dlab <- data.frame(
    I =  c(0.0,0.1,0.2,0.3,0.4,0.6,0.8),
    xB = c(0.0,0.8,0.8,0.8,0.8,0.85,0.9)
  )
  for (i in 1:nrow(dlab)){
    dlab$xB[i] <- unique(dl$xB)[which.min(abs(dl$xB-dlab$xB[i]))]  
    dlab$yB[i] <- approx(x=rev(dl$I[dl$xB==dlab$xB[i]]), y=rev(dl$zB[dl$xB==dlab$xB[i]]), xout=dlab$I[i])$y
  }
  dlab$xB[dlab$I==0] <- 0.5*(min(xlim))
  dlab$yB[dlab$I==0] <- 0.0
  dlab$label <- dlab$I
  dlab$label[dlab$I==0.0] <- 'I=0.0'

  #plot
  p <- ggplot() + 
    theme_bw() + 
    theme(legend.position='none') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  if (showgrid==T){
    p <- p + geom_tile(data=dl, aes(x=xB, y=zB, fill=I))
  }
  p <- p + 
    geom_hline(yintercept=seq(-5,10,by=0.5), size=0.2, color='grey90') +
    geom_vline(xintercept=seq(-5,5,by=0.5), size=0.2, color='grey90') +
    geom_hline(yintercept=0) + 
    stat_contour(data=dl, aes(x=xB, y=zB, z=I), binwidth=0.1, color='black', linetype=1) + 
    xlab("x/B [-]") +
    ylab("z/B [-]") +
    scale_y_reverse() + #
    coord_fixed(ratio=1, xlim=xlim, ylim=c(ylim[2],ylim[1]-2*larrow), expand=F) + 
    scale_fill_distiller(palette='Blues', direction=1) +
    annotate('segment',x=0.2,y=-0.2*larrow,xend=0.2,yend=0, arrow=arrow(length=unit(0.03,"npc"))) +
    annotate('segment',x=0.4,y=-0.4*larrow,xend=0.4,yend=0, arrow=arrow(length=unit(0.03,"npc"))) +
    annotate('segment',x=0.6,y=-0.6*larrow,xend=0.6,yend=0, arrow=arrow(length=unit(0.03,"npc"))) +
    annotate('segment',x=0.8,y=-0.8*larrow,xend=0.8,yend=0, arrow=arrow(length=unit(0.03,"npc"))) +
    annotate('segment',x=1.0,y=-1.0*larrow,xend=1.0,yend=0, arrow=arrow(length=unit(0.03,"npc"))) +
    annotate('segment',x=0.0,y=0.0,xend=1.0,yend=-larrow) +
    geom_label(data=dlab, aes(x=xB,y=yB,label=label), fill='white', alpha=0.7, size=labelsize,
               label.padding = unit(0.1, "lines"),
               label.r = unit(0.1, "lines"))
  
  #return
  return(p)
}


#add points on chart
f_trianglestresschart_adddata <- function(
  p,
  d,
  legpos = NA,
  crosshairs = T,
  values = T,
  nround = 2
)
{
  #copy chart
  p2 <- p
  #add labels
  if (!'label'%in%colnames(d)){
    d$label <- as.factor(seq(nrow(d)))
  }
  #calculate I
  fI <- function(xB,zB){
    delta <- atan2(xB-1,zB)
    alphadelta <- atan2(xB,zB)
    alpha <- alphadelta - delta
    I <- (xB*alpha - 0.5*sin(2*delta))/pi
    return(I)
  }
  d$I <- fI(d$xB, d$zB)
  #add crosshairs
  if (crosshairs==T){
    dl <- d %>%
      group_by(label) %>%
      group_modify(~data.frame(xB=c(-10,.x$xB,.x$xB),zB=c(.x$zB,.x$zB,10)))
    p2 <- p + 
      geom_path(data=dl, aes(x=xB,y=zB,color=as.factor(label)))
  }
  #add values
  if (values==T){
    d$Ilab <- with(d, paste("I==",round(I,nround),sep=''))
    p2 <- p2 + 
      geom_text(data=d, aes(x=xB+0.1, y=zB, label=Ilab, color=as.factor(label)), hjust=0, vjust=0.5, parse=T, show.legend=F)
  }
  #add data points 
  p2 <- p2 + 
    geom_point(data=d, aes(x=xB, y=zB, shape=as.factor(label), color=as.factor(label))) +
    scale_color_brewer(name='',palette='Set1')
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

