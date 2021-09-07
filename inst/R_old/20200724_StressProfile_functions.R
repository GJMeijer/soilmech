#Load packages
library('plotly')
library('tidyverse')
source('rscripts/R_StandardFunctions.R')

#functions
# - ds <- f_stress(d)
# - f_plotly_together(d)
#   - f_plotly_soilprofile(d)
#   - f_plotly_stressprofile(ds)
# - f_ggplot_stressprofile(d)


#function to calculate stresses + labels
f_stress <- function(
  d,     #dataframe with data:
         # - z: depth at which new layer starts
         # - q: surcharge applied at this point 
         # - saturated: layer saturated with water (true of false)
         # - gammab:  bulk density
  zmax=10,   #maximum depth
  marg=1e-6,   #margin parameters, does not do much else
  zinterval=NA,   #add points at regular interval. parameter specifies distance
  labelmode='stressonly',   #plotly labels
  gammaw=10,  #unit weigth of water
  output='plot'   #if 'bookdown', generate neat bookdown table and header names
){
  #add points on regular interval if required
  if (!is.na(zinterval)) {
    #create points with default values
    dadd <- data.frame(
      z=seq(seq(min(d$z), max(d$z), zinterval)),
      q = 0,
      saturated = T,
      gammab = 20
    )
    #remove points already defined in <d>
    dadd <- dadd[!dadd$z%in%d$z,]
    #change saturation and unit weight when required
    for (i in 1:dim(dadd)[1]){
      ind <- min(which(d$z>=dadd$z[i])[1]-1,dim(d)[1],na.rm=T)
      dadd$saturated[i] <- d$saturated[ind]
      dadd$gammab[i] <- d$gammab[ind]
    }
    #merge and sort
    d <- rbind(d,dadd)
    d <- d[order(d$z),]
  }
  #get dataframe with stresses - excluding points just before surcharge
  ds1 <- data.frame(
    layer = 1+seq(0,dim(d)[1]),
    z = c(d$z,zmax),
    sigma = cumsum(c(0,diff(c(d$z,zmax))*d$gammab)) + cumsum(c(d$q,0)),
    u = cumsum(c(0,d$saturated*diff(c(d$z,zmax))*gammaw))
  )
  #increase with respect to last point
  ds1$dgammab <- c(0,d$gammab)
  ds1$dq <- c(d$q,0)
  ds1$dz <- c(0,diff(c(d$z,zmax)))
  ds1$sigma_last <- c(0,head(ds1$sigma,-1))
  ds1$u_last <- c(0,head(ds1$u,-1))
  #get dataframe with stresses - just before surcharge
  ds0 <- ds1[c(d$q!=0,F),]
  ds0$z <- ds0$z - marg
  ds0$sigma <- ds0$sigma - d$q[d$q!=0]
  #increase with espect to last point
  if (dim(ds0)[1]>0){
    ds0$dq <- 0
  }
  #merge together and sort
  ds <- rbind(ds1, ds0[ds0$z>=min(ds1$z),])
  ds <- ds[order(ds$z),]
  #effective stress
  ds$sigmad <- with(ds, sigma-u)

  #labels
  ds$sigma_label <- NA
  ds$u_label <- NA
  ds$sigmad_label <- NA
  if (labelmode=='stressonly'){
    ds$sigma_label <- with(ds, paste('\u03c3<sub>v</sub> = ', ds$sigma, ' kPa'))
    ds$u_label <- with(ds, paste('u = ', ds$sigma, ' kPa'))
    ds$sigmad_label <- with(ds, paste("\u03c3'<sub>v</sub> = ", ds$sigmad, ' kPa'))
  } else if (labelmode=='calculations'){
    #labels - first point
    if (ds$dq[1]==0) {
      ds$sigma_label[1] <- '\u03c3<sub>v</sub> = 0 kPa'
    } else {
      ds$sigma_label[1] <- paste('\u03c3<sub>v</sub> = q = ', ds$dq[1], ' kPa', sep='')
    }
    ds$u_label[1] <- 'u = 0 kPa'
    ds$sigmad_label[1] <- paste("\u03c3'<sub>v</sub> = \u03c3<sub>v</sub> - u = ", ds$sigma[1], ' - ', ds$u[1], ' = ', ds$sigmad[1], ' kPa', sep='')
    # other labels - total
    ind0 <- which(ds$dq==0 & ds$layer>1)
    ds$sigma_label[ind0] <- with(ds[ind0,], paste('\u03c3<sub>v</sub> = ', sigma_last , ' + ', dz, '\u00b7', dgammab, ' = ', sigma, ' kPa', sep=''))
    ind1 <- which(ds$dq!=0 & ds$layer>1)
    #ds$sigma_label[ind1] <- with(ds[ind1,], paste('\u03c3 = ', sigma_last , ' + ', dz, '\u00b7', dgammab, ' + ', dq, ' = ', sigma, ' kPa', sep=''))
    ds$sigma_label[ind1] <- with(ds[ind1,], paste('\u03c3<sub>v</sub> = ', ds$sigma[ind1-1], ' + ', dq, ' = ', sigma, ' kPa', sep=''))
    # labels - pore water pressure
    ds$u_label[ds$u==0] <- 'u = 0 kPa'
    ds$u_label[ds$u>0] <- with(ds[ds$u>0,], paste('u = ', u_last, ' + ', dz, '\u00b710 = ', u, ' kPa', sep=''))
    # labels - effective stress
    ds$sigmad_label <- with(ds, paste("\u03c3'<sub>v</sub> = \u03c3<sub>v</sub> - u = ", sigmad, ' kPa', sep=''))
  }
  #return
  if (output=='bookdown'){
    return(list(
      ds[,c('z','sigma','u','sigmad')], 
      c('Depth [m]','$\\sigma_v$ [kPa]','$u$ [kPa]',"$\\sigma'_v$ [kPa]")
    ))
  } else {
    return(ds)
  }
}


#plot function for stress profiles
f_plotly_stressprofile <- function(
  ds,
  linewidth=1.5,
  cols=c('#000000','#0000ff','#ff0000'),
  linetype=c('solid','dash','dot'),
  xlim=NA,
  ylim=NA,
  traces=c('sigma','u','sigmad'),
  hovermode='y'
){

  #axis limits
  if (is.na(xlim[1])){
    xlim = c(0,max(ds$sigma))
  } else {
    xlim <- sort(xlim,'decreasing'=T)
  }
  if (is.na(ylim[1])){
    ylim = c(max(ds$z), min(ds$z))
  } else {
    ylim <- sort(ylim,'decreasing'=T)
  }

  #initiate plot
  p <- plot_ly(ds)
  #add traces
  if ('sigma'%in%traces){
    p <- p %>% add_trace(
      x = ~sigma,
      y = ~z,
      name = 'Total stress',
      type='scatter',
      mode = 'lines+markers',
      hoverinfo = 'text',
      text = ~sigma_label,
      marker = list(
        color=cols[1]
      ),
      line = list(
        width=linewidth,
        color=cols[1],
        dash=linetype[1]
      )
    )
  }
  if ('u'%in%traces){
    p <- p %>% add_trace(
      x = ~u,
      y = ~z,
      name = 'Pore water pressure',
      type='scatter',
      mode = 'lines+markers',
      hoverinfo = 'text',
      text = ~u_label,
      marker = list(
        color=cols[2]
      ),
      line = list(
        width=linewidth,
        color=cols[2],
        dash=linetype[2]
      )
    )
  }
  if ('sigmad'%in%traces){
    p <- p %>% add_trace(
      x = ~sigmad,
      y = ~z,
      name = 'Effective stress',
      type='scatter',
      mode = 'lines+markers',
      hoverinfo = 'text',
      text = ~sigmad_label,
      marker = list(
        color=cols[3]
      ),
      line = list(
        width=linewidth,
        color=cols[3],
        dash=linetype[3]
      )
    )
  }
  #add layout
  p <- p %>% layout(
    xaxis = list(
      range=xlim,
      title='Stress [kPa]',
      side='top',
      visible=T
    ),
    yaxis = list(
      range=ylim,
      title='Depth [m]'
    ),
    legend = list(
      x=1,
      y=1,
      xanchor='right',
      yanchor='top'
    ),
    hovermode = hovermode,
    template = 'none'
  )
  #changes axis titles if only one plot
  if (length(traces)==1){
    if (traces=='sigma'){
      p <- p %>% layout(xaxis=list(title='Total stress [kPa]'))
    } else if (traces=='u'){
      p <- p %>% layout(xaxis=list(title='Pore water pressure [kPa]'))
    } else if (traces=='sigmad'){
      p <- p %>% layout(xaxis=list(title='Effective stress [kPa]'))
    }
  }
  #return
  return(p)
}

#plotly function for soil profiles
f_plotly_soilprofile <- function(
  d,
  zmax = 10,
  gammaw = 10,
  col_soildry = '#d3bc5f',
  col_soilsat = '#aebab7',
  col_soilline = '#65571d',
  col_water = '#2a7fff',
  col_waterline = '#2a7fff',
  linetype_soil = 'solid',
  linetype_water = 'dash',
  linewidth = 2.0,
  xlim = NA,
  ylim = NA,
  arrow_size = 0.075,
  arrow_n = 6){

  #axis limits
  if (is.na(ylim[1])){
    ylim <- c(max(d$z),min(d$z))
  }
  if (is.na(xlim[1])){
    xlim <- c(0,1)
  }

  #initiate plot
  p <- plot_ly()
  #add layers + annotation
  for (i in 1:(dim(d)[1])) {
    #get coordinates
    if (i==dim(d)[1]) {
      z0 <- d$z[i]
      z1 <- zmax
    } else {
      z0 <- d$z[i]
      z1 <- d$z[i+1]
    }
    #hover labels for soil masses
    label1 <- paste('\u03B3<sub>b</sub> = ', d$gammab[i], ' kN/m\u00b3', sep='')
    label2 <- paste('h = ', c(d$z[-1],zmax)-d$z, ' m', sep='')
    label3 <- rep('Fully saturated soil', dim(d)[1])
    label3[d$saturated==F] <- 'Dry soil'
    label3[d$saturated==T & d$gammab==gammaw] <- 'Water'
    label <- paste(label3,label1,label2,sep='<br>')
    #get color
    if (d$saturated[i]==F){
      soil_color <- col_soildry
    } else if (d$saturated[i]==T & d$gammab[i]==gammaw) {
      soil_color <- col_water
    } else {
      soil_color <- col_soilsat
    }
    #plot rectangle
    p <- p %>% add_trace(
      type = 'scatter',
      mode = 'lines',
      x = c(xlim, rev(xlim)),
      y = c(z0,z0,z1,z1),
      fill = 'toself',
      fillcolor = soil_color,
      line = list('width'=0.0),
      hoveron = 'fills',
      text = label[i], #paste('\u03B3<sub>d</sub> = ', d$gammab[i], ' kN/m\u00b3<br>h = ', z1-z0, ' m', sep=''),
      hoverinfo='text',
      showlegend=F
    )
    #add text
    p <- p %>% add_trace(
      type='scatter',
      x=0.5*(xlim[2]-xlim[1]),
      y=0.5*(z0+z1),
      mode="text",
      text=paste('\u03B3<sub>b</sub> = ', d$gammab[i], ' kN/m\u00b3', sep=''),
      textposition="middle center",
      showlegend=F,
      hoverinfo='skip'
    )
    #add interface
    if (c(d$saturated[1], diff(d$saturated))[i]==1){
      p <- p %>% add_trace(
        type='scatter',
        mode = 'lines',
        x = xlim,
        y = c(z0,z0),
        line = list(
          width=linewidth,
          color=col_waterline,
          dash=linetype_water
        ),
        hoverinfo='skip',
        showlegend=F
      )
    } else {
      p <- p %>% add_trace(
        type='scatter',
        mode = 'lines',
        x = xlim,
        y = c(z0,z0),
        line = list(
          width=linewidth,
          color=col_soilline,
          dash=linetype_soil
        ),
        hoverinfo='skip',
        showlegend=F
      )
    }
    #add overburden arrows
    if (d$q[i]!=0){
      #arrow data
      da <- data.frame(
        x0=seq(xlim[1],xlim[2]-1/arrow_n,l=arrow_n) + 0.5/arrow_n,
        x1=seq(xlim[1],xlim[2]-1/arrow_n,l=arrow_n) + 0.5/arrow_n,
        y0=z0,
        y1=z0-(max(d$z)-min(d$z))*arrow_size
      )
      #add to plot
      p <- p %>% add_annotations(
        x = da$x0,
        y = da$y0,
        xref = "x", yref = "y",
        axref = "x", ayref = "y",
        text = "",
        showarrow = T,
        ax = da$x1,
        ay = da$y1
      )
      #add overburden text
      p <- p %>% add_trace(
        type= 'scatter',
        mode='text',
        x = xlim[1], #0.5*(xlim[1]+xlim[2]),
        y = da$y1[1],
        text = paste('q = ', d$q[i], ' kPa', sep=''),
        textposition="top right",
        hoverinfo='skip',
        showlegend=F
      )
      #change axis limits
      ylim[2] <- min(ylim[2], z0-2*(max(d$z)-min(d$z))*arrow_size)
    }
  }
  #add layout
  p <- p %>% layout(
    xaxis = list(
      range = xlim,
      title = 'position',
      side = 'top',
      visible = F
    ),
    yaxis = list(
      range = ylim,
      title = 'Depth [m]'
    ),
    showlegend = F,
    template = 'none'
  )
  #return
  return(p)
}

#plot together
f_plotly_together <- function(d,
                              zmax=10,
                              traces=c('sigma','u','sigmad'),
                              zinterval=NA,
                              labelmode='stressonly'){
  ds <- f_stress(d, zmax=zmax, zinterval=zinterval, labelmode=labelmode)
  p1 <- f_plotly_soilprofile(d, zmax=zmax)
  p2 <- f_plotly_stressprofile(ds, traces=traces)
  p <- subplot(p1,p2,
               shareY=T, titleX=T,
               widths=c(0.3,0.7),
               margin=0.05)
  p$width <- 600
  p$height <- 300
  return(p)
}


#ggplot profile
f_ggplot_stressprofile <- function(
  ds,
  linewidth=0.5,
  cols=c('#000000','#0000ff','#ff0000'),
  linetype=c(1,2,3),
  xlim=NA,
  ylim=NA,
  zmax=10,
  traces=c('sigma','u','sigmad'),
  layers=T
){
  #calculate stresses
  ds <- f_stress(d, zmax=zmax)
  #initiate plot
  p <- ggplot()
  #convert to long data
  dsl <- gather(ds[,c('z',traces)],stresstype,value,traces)
  
  #add data
  p <- p + geom_path(data=dsl, aes(y=z, x=value, color=stresstype, linetype=stresstype))
  if (layers==T){
    p <- p + geom_hline(data=ds, aes(yintercept=z), color='grey50', linetype=4)
  }
  #select axis labels
  if (length(traces)==1){
    if (traces=='sigma'){
      xlab <- 'Total stress [kPa]'
    } else if (traces=='u') {
      xlab <- 'Pore water pressure [kPa]'
    } else if (traces=='sigmad') {
      xlab <- 'Effective stress [kPa]'
    } else {
      xlab <- 'Stress [kPa]'  #technicall not required, but just in case of spelling error or something like that
    }
    legend <- F
  } else {
    xlab <- 'Stress [kPa]'
    legend <- T
  }
  #axis extend
  if (is.na(xlim)) {
    xlim <- c(0,f_roundlimits(max(dsl$value)))
  }
  if (is.na(ylim)) {
    ylim <- c(max(ds$z), min(ds$z))
  }
  #change axes
  p <- p + 
    theme_bw() +
    xlab(xlab) +
    ylab('Depth [m]') +
    scale_x_continuous(lim=xlim, expand=c(0,0), position='top') + 
    scale_y_reverse(lim=ylim, expand=c(0,0))
  #selected traces
  sel <- sort(match(traces,c('sigma','u','sigmad')))
  #change colors
  p <- p + 
    scale_linetype_manual(name='', values=linetype[sel], labels=expression(sigma[v],u,sigma*"'"[v])[sel]) + 
    scale_color_manual(name='', values=cols[sel], labels=expression(sigma[v],u,sigma*"'"[v])[sel])
  #legend
  if (legend==T){
    p <- p + theme(
      legend.position=c(0.95,0.95), 
      legend.justification=c(1,1),
      legend.title=element_blank()
    )
  } else {
    p <- p + theme(legend.position='none')
  }
  #return
  return(p)
}


 #test
#d <- data.frame(
#  z = c(0,5,6),
#  q = c(50,0,0),
#  gammab = c(20,10,18),
#  saturated = c(F,T,T)
#)
#ds <- f_stress(d, output='bookdown')
#f_ggplot_stressprofile(d, traces=c('sigmad','sigma'))

