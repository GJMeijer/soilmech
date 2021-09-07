library('tidyverse')
library("RColorBrewer")
library('directlabels')

#empty plot
f_stripstresschart_empty <- function(
  xlim = c(-2.5,2.5),
  ylim = c(0,6.5),
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
    delta <- atan2(xB-0.5,zB)
    alphadelta <- atan2(xB+0.5,zB)
    alpha <- alphadelta - delta
    I <- (alpha+sin(alpha)*cos(alpha+2*delta))/pi
    return(I)
  }
  
  #cals
  dl$I <- fI(dl$xB,dl$zB)

  #labels
  dlab <- data.frame(
    I = c(0,0.1,0.2,0.3,0.4,0.6,0.8,1.0),
    xB = 0
  )
  dlab$yB = approx(x=rev(dl$I[dl$xB==0]), y=rev(dl$zB[dl$xB==0]), xout=dlab$I)$y
  dlab$xB[dlab$I==0] <- 0.5*(min(xlim)-0.5)
  dlab$yB[dlab$I==0] <- 0.0
  dlab$label <- dlab$I
  dlab$label[dlab$I==0.0] <- 'I=0.0'

  #plot
  p <- ggplot() + 
    theme_bw() + 
    theme(legend.position='none') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  if (showgrid==T){
    p <- p + geom_tile(data=dl, aes(x=xB, y=zB, fill=I), show.legend=F)
  }
  p <- p + 
    geom_hline(yintercept=seq(-5,10,by=0.5), size=0.2, color='grey90') +
    geom_vline(xintercept=seq(-5,5,by=0.5), size=0.2, color='grey90') +
    geom_hline(yintercept=0) + 
    stat_contour(data=dl, aes(x=xB, y=zB, z=I), binwidth=0.1, color='black', linetype=1, show.legend=F) + 
    xlab("x/B [-]") +
    ylab("z/B [-]") +
    scale_y_reverse() + #
    coord_fixed(ratio=1, xlim=xlim, ylim=c(ylim[2],ylim[1]-2*larrow), expand=F) + 
    scale_fill_distiller(palette='Blues', direction=1) +
    annotate('segment',x=-0.5,y=-larrow,xend=-0.5,yend=0, arrow=arrow(length=unit(0.03,"npc"))) +
    annotate('segment',x=-0.3,y=-larrow,xend=-0.3,yend=0, arrow=arrow(length=unit(0.03,"npc"))) +
    annotate('segment',x=-0.1,y=-larrow,xend=-0.1,yend=0, arrow=arrow(length=unit(0.03,"npc"))) +
    annotate('segment',x=0.1,y=-larrow,xend=0.1,yend=0, arrow=arrow(length=unit(0.03,"npc"))) +
    annotate('segment',x=0.3,y=-larrow,xend=0.3,yend=0, arrow=arrow(length=unit(0.03,"npc"))) +
    annotate('segment',x=0.5,y=-larrow,xend=0.5,yend=0, arrow=arrow(length=unit(0.03,"npc"))) +
    annotate('segment',x=-0.5,y=-larrow,xend=0.5,yend=-larrow) +
    geom_label(data=dlab, aes(x=xB,y=yB,label=label), fill='white', alpha=0.7, size=labelsize,
               label.padding = unit(0.1, "lines"),
               label.r = unit(0.1, "lines"))
  
  #return
  return(p)
}

#add points on chart
f_stripstresschart_adddata <- function(
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
    delta <- atan2(xB-0.5,zB)
    alphadelta <- atan2(xB+0.5,zB)
    alpha <- alphadelta - delta
    I <- (alpha+sin(alpha)*cos(alpha+2*delta))/pi
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
