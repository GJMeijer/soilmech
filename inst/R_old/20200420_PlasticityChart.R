#Load packages
library('tidyverse')
#library('RColorBrewer')


#functions
# - f_create_empty_plasticity_chart()
# - f_create_plasticity_chart(d)  
#   - dataframe <d> must have fields <wL>,<Ip> and optionally <label>
# - f_createplasticityfit(d)
#   - dataframe <d> must have fields <w> and <u>


#############
### INPUT ###
#############

#Function to create empty plasticity chart
f_create_empty_plasticity_chart <- function(label=T){
  # INPUT
  # - <label>: if True, adds legend on top left for symbols in plot
  
  #Plasticity limits
  dlim <- data.frame(Symbol=c('L','I','H','V','E'),
                     Text=c('Low plasticity','Medium plasticity','High plasticity','Very high plasticity','Extremely high plasticity'))
  wLmin <- 0
  wLmax <- 120
  Ipmax <- 80
  dlim$wLmin <- c(wLmin,35,50,70,90)
  dlim$wLmax <- c(35,50,70,90,wLmax) 
  dlim$wLavg <- with(dlim, 0.5*(wLmin+wLmax))
  
  #A-B-line
  dAB <- data.frame(wL=c(20,20+6.3/0.73,wLmax, 15,wLmax))
  dAB$Type <- c('A','A','A','B','B')
  dAB$Ip <- pmax(6.3, 0.73*(dAB$wL-20))
  dAB$Ip[dAB$Type=='B'] <- with(dAB[dAB$Type=='B',], pmax(6.3, 0.9*(wL-8)))
  
  #Symbols
  dsym <- data.frame(PSymbol=rep(dlim$Symbol,2))
  dsym$SSymbol <- c(rep('M',5),rep('C',5))
  dsym$Symbol <- with(dsym,paste(SSymbol,PSymbol,sep=''))
  dsym$wL <- pmax(30,rep(dlim$wLavg,2))
  dsym$IpA <- pmax(6.3, 0.73*(dsym$wL-20))
  dsym$IpB <-  pmax(6.3, 0.9*(dsym$wL-8))
  dsym$Ip <- 0.5*dsym$IpA
  dsym$Ip[dsym$SSymbol=='C'] <- with(dsym[dsym$SSymbol=='C',],0.5*(IpA+pmin(Ipmax,IpB)))
  
  #Plot settings
  textsize1 <- 3
  textsize2 <- 2.5
  lab <- 'Predominant behaviour\nmaterial <0.425 mm\n- C = Clay\n- M = Silt\n\nPlasticity:\n- L = Low\n- I = Intermediate\n- H = High\n- V = Very High\n- E = Extremely high'
  #Create plot
  p <- ggplot() + 
    theme_bw() + theme(legend.position='none') + 
    geom_line(data=dAB, aes(x=wL, y=Ip, linetype=Type), show.legend=F) +
    geom_vline(data=dlim, aes(xintercept=wLmin), size=0.5, linetype=3, show.legend=F) +
    scale_x_continuous(breaks=seq(0,wLmax,20)) +
    scale_y_continuous(breaks=seq(0,Ipmax,20)) +
    coord_cartesian(xlim=c(wLmin,wLmax), ylim=c(0,Ipmax), expand=F) + 
    xlab(expression(Liquid ~limit~w[L]~"[%]")) + 
    ylab(expression(Plasticity~index~I[p]~"[%]")) +
    geom_text(data=dsym, aes(x=wL,y=Ip,label=Symbol), show.legend=F) + 
    annotate('label', x=105, y=0.73*(105-20), label="'A'-line", size=textsize1) + 
    annotate('label', x=15, y=0.90*(15-8), label="'B'-line", hjust=1, size=textsize1) 
  if (label==T){
    p <- p + annotate('label', x=wLmin+0.01*(wLmax-wLmin),y=0.99*Ipmax, label=lab, size=textsize2, hjust=0, vjust=1)
  }
  #return
  return(p)
}

#function to add points to plasticity chart
f_create_plasticity_chart <- function(d, legpos=NA, crosshairs=T, label=F){
  # INPUT
  # - <d>: dataframe with fields <wL>, <Ip> and legend labels <label> if requested
  # - <legpos>: legend position. If NA, do not use legend
  # - <crosshairs> if true, add lines pointing towards points
  #create empty plot
  p <- f_create_empty_plasticity_chart(label=label)
  #if no label
  if (!'label'%in%colnames(d)){
    d$label <- as.factor(seq(nrow(d)))
  }
  #add crosshairs
  if (crosshairs==T){
    dl <- d %>%
      group_by(label) %>%
      group_modify(~data.frame(wL=c(.x$wL,.x$wL,0),Ip=c(0,.x$Ip,.x$Ip)))
    p <- p + 
      geom_path(data=dl, aes(x=wL,y=Ip,color=as.factor(label), linetype=as.factor(label)))
  }
  #add points
  p <- p + 
    geom_point(data=d, aes(x=wL, y=Ip, color=as.factor(label), shape=as.factor(label)), size=2) + 
    scale_color_brewer(name='Sample', palette='Set1') 
  #legend
  if (is.na(legpos)){
    p <- p + theme(legend.position='none')
  } else {
    p <- p + theme(legend.position=legpos)
  }
  #return
  return(p)
}


#function to create linear fit through data - plasticity
f_createplasticityfit <- function(d=NA, xlim=c(40,80), ylim=c(0,30), 
                                  arrowlength=0.2, label_size=4, nsignif=2, palette='Set1'){
  # INPUT
  # - <d>: dataframe with fields <u> and <w>
  #        if not given, plot empty chart
  #fit data
  if (!is.na(d)){
    ft <- lm(u~w, data=d)
    wL <- (20-coef(ft)[1])/coef(ft)[2]
    #create plot
    p <- ggplot(d, aes(x=w,y=u)) +
      theme_bw() +
      geom_smooth(method='lm', se=F, size=0.5, linetype=2, color='darkgreen') +
      geom_point(color='darkblue', size=2) + 
      coord_cartesian(xlim=xlim, ylim=ylim, expand=F) + 
      annotate('segment', x=0, y=20, xend=wL, yend=20, arrow=arrow(length=unit(arrowlength,"cm")), color='red') + 
      annotate('segment', x=wL, y=20, xend=wL, yend=0, arrow=arrow(length=unit(arrowlength,"cm")), color='red') +
      scale_color_brewer(palette=palette) + 
      annotate('text', x=wL, y=20, hjust=0, vjust=0.5, label=paste('~w[L]==',signif(wL,nsignif),'*"%"',sep=''), parse=T, size=label_size) +
      xlab('w [%]') +
      ylab('penetration depth [mm]')
  } else {
    p <- ggplot() +
      theme_bw() +
      coord_cartesian(xlim=xlim, ylim=ylim, expand=F) + 
      xlab('w [%]') +
      ylab('penetration depth [mm]')
  }
  #return plot and value
  return(list(signif(wL,nsignif), p))
}
