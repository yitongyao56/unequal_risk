load(file='/Users/yyao/Downloads/Cscat/0127radar_global_grid_hydro_annual_1992-2022.Rdata')
global.radar.hydro.annual.mean=apply(global.radar.hydro.annual[, , (1993-1993+1):(2020-1993+1)], c(1,2), mean, na.rm=TRUE) 

library(raster)
projlat<-CRS("+proj=longlat +datum=WGS84 +no_defs")
amazon.Rmean<-raster(nrow=600,ncol=600)  # 0.05 deg 
extent(amazon.Rmean)<-c(-80, -50,-20,10)
amazon.Rmean[is.na(amazon.Rmean)]<-1
projection(amazon.Rmean)<-projlat
regrid_cscat_chirps_amazon=function(mat) {
    #reproject1d<-function(mat,Rmean){
    pr_mask<-raster(t(mat[, 120:1]))
    extent(pr_mask)<-c(-80,-50,-20,10)
    projection(pr_mask)<-projlat
        
    meanpr_proj<-projectRaster(from = pr_mask,to=amazon.Rmean)
    mat_pr<-t(as.matrix(meanpr_proj))[,600:1]
    return(mat_pr)
    #}
}


congo.Rmean<-raster(nrow=400,ncol=600)  # 0.05 deg 
extent(congo.Rmean)<-c(5,35,-10,10)
congo.Rmean[is.na(congo.Rmean)]<-1
projection(congo.Rmean)<-projlat
regrid_cscat_chirps_congo=function(mat) {
    #reproject1d<-function(mat,Rmean){
    pr_mask<-raster(t(mat[, 80:1]))
    extent(pr_mask)<-c(5,35,-10,10)
    projection(pr_mask)<-projlat
        
    meanpr_proj<-projectRaster(from = pr_mask,to=congo.Rmean)
    mat_pr<-t(as.matrix(meanpr_proj))[,400:1]
    return(mat_pr)
    #}
}

tiff('elnino_basin_cscat_anomaly Oct17.tiff', width=850, height=800)
par(mfrow=c(2,2))
label.num=c('(a)','(b)')
label.year=c('1997/98','2015/16')
xx=0
for (year in c(1998, 2016)) {
    global.radar.hydro.annual.anomaly=global.radar.hydro.annual[, , year-1993+1]-global.radar.hydro.annual.mean

    amazon.radar.anomaly=global.radar.hydro.annual.anomaly[401:520, 281:400]  # -20 
    #africa.radar.anomaly=global.radar.hydro.annual.anomaly[761:860, 321:400]  # 10 190*4=760 

    draw.anomaly=regrid_cscat_chirps_amazon(amazon.radar.anomaly)
    draw.anomaly[amazon.tmf.005<0.95]=NA
    draw.anomaly[draw.anomaly>0.25]=0.25
    draw.anomaly[draw.anomaly<(-0.25)]=-0.25
    colos=brewer.pal(10,'RdYlGn')

    pie.amazon=array(NA, 10)
    bk.bins=seq(-0.25, 0.25, by=0.05)
    bk.bins[11]=0.26
    for (bb in 1:10) {
        temp=((draw.anomaly>=bk.bins[bb]) & (draw.anomaly<bk.bins[bb+1]))*amazon.area
        pie.amazon[bb]=sum(temp, na.rm=TRUE)
    }

    pie.amazon.frac=pie.amazon/sum((!is.na(draw.anomaly))*amazon.area, na.rm=TRUE)

    ##
    if (year==1998) {
        par(fig=c(0.05, 0.5, 0.52, 0.975), mar=c(0,0,0,0), new=TRUE)  
    } else {
        par(fig=c(0.525, 0.975, 0.52, 0.975), mar=c(0,0,0,0), new=TRUE) 
    }
    par(mar=c(5,3,3,3))
    xx=xx+1
    image(seq(-80,-50,by=1/20), seq(-20, 10, by=1/20), draw.anomaly, xlab='', ylab='', cex.axis=2, cex.lab=2, las=1,
    axis.args=list(cex.axis=2, cex=2), zlim=c(-0.25, 0.25), col=colos, legend.shrink=1, main=paste0(label.num[xx], ' Amazon ', label.year[xx]), cex.main=2, 
    breaks=seq(-0.25, 0.25, by=0.05), lab.breaks=c('≤-0.25','-0.2','-0.15','-0.10','-0.05','0','0.05','0.10','0.15','0.20', '≥0.25'))
    map(database = "world",
            mar=c(1,1,1,1),col='gray', panel.first = grid(),add=TRUE,lwd=2)
    map.axes(cex.axis=2,lwd=2, las=1)

    if (year==1998) {
        par(fig=c(0.1, 0.2, 0.58, 0.7), mar=c(0,0,0,0), new=TRUE)  
    } else {
        par(fig=c(0.575, 0.675, 0.58, 0.7), mar=c(0,0,0,0), new=TRUE) 
    }
    pie(pie.amazon.frac,,col=colos, labels=NA)
}
label.num=c('(c)','(d)')
xx=0
for (year in c(1998,2016)) {
    global.radar.hydro.annual.anomaly=global.radar.hydro.annual[, , year-1993+1]-global.radar.hydro.annual.mean

    #amazon.radar.anomaly=global.radar.hydro.annual.anomaly[401:520, 281:400]  # -20 
    africa.radar.anomaly=global.radar.hydro.annual.anomaly[741:860, 321:400]  # 10 190*4=760 start from -180 ends 5 deg 

    draw.anomaly=regrid_cscat_chirps_congo(africa.radar.anomaly)
    draw.anomaly[africa.tmf.005.match<0.95]=NA
    draw.anomaly[draw.anomaly>0.25]=0.25
    draw.anomaly[draw.anomaly<(-0.25)]=-0.25

    pie.africa=array(NA, 10)
    bk.bins=seq(-0.25, 0.25, by=0.05)
    bk.bins[11]=0.26
    for (bb in 1:10) {
        temp=((draw.anomaly>=bk.bins[bb]) & (draw.anomaly<bk.bins[bb+1]))*congo.area
        pie.africa[bb]=sum(temp, na.rm=TRUE)
    }

    pie.africa.frac=pie.africa/sum((!is.na(draw.anomaly))*congo.area, na.rm=TRUE)

    colos=brewer.pal(10,'RdYlGn')
    if (year==1998) {
        par(fig=c(0.05, 0.5, 0.07, 0.52), mar=c(0,0,0,0), new=TRUE)  
    } else {
        par(fig=c(0.525, 0.975, 0.07, 0.52), mar=c(0,0,0,0), new=TRUE) 
    }
    par(mar=c(5,3,3,3))
    xx=xx+1
    image(seq(5,35,by=1/20), seq(-10, 10, by=1/20), draw.anomaly, xlab='', ylab='', cex.axis=2, cex.lab=2, xlim=c(7,35),las=1,
    axis.args=list(cex.axis=2, cex=2), zlim=c(-0.25, 0.25), col=colos, legend.shrink=1, main=paste0(label.num[xx], ' C. Africa ',label.year[xx]), cex.main=2, 
    breaks=seq(-0.25, 0.25, by=0.05), lab.breaks=c('≤-0.25','-0.2','-0.15','-0.10','-0.05','0','0.05','0.10','0.15','0.20', '≥0.25'))
    map(database = "world",
            mar=c(1,1,1,1),col='gray', panel.first = grid(),add=TRUE,lwd=2)
    map.axes(cex.axis=2,lwd=2, las=1)

    if (year==1998) {
        par(fig=c(0.1, 0.2, 0.13, 0.25), mar=c(0,0,0,0), new=TRUE)  
    } else {
        par(fig=c(0.575, 0.675, 0.13, 0.25), mar=c(0,0,0,0), new=TRUE) 
    }
    pie(pie.africa.frac,,col=colos, labels=NA)
}

par(fig=c(0.15, 0.87, 0.08, 0.25), mar=c(0,0,0,0), new=TRUE)
image.plot(legend.only=TRUE, horizontal=TRUE, zlim=c(-0.25, 0.25), legend.lab=expression(Delta~'BK (dB)'), legend.cex=2,  axis.args=list(cex.axis=2, cex=2), 
col=colos, breaks=seq(-0.25, 0.25, by=0.05), legend.shrink=1, legend.width=3, legend.line=3, 
lab.breaks=c('≤-0.25','-0.2','-0.15','-0.10','-0.05','0','0.05','0.10','0.15','0.20', '≥0.25'))

dev.off()
