
load(file='/Users/yyao/Downloads/Cscat/0127radar_global_grid_hydro_annual_1992-2022.Rdata')
global.radar.hydro.annual.mean=apply(global.radar.hydro.annual[, , (1993-1993+1):(2020-1993+1)], c(1,2), mean, na.rm=TRUE) 

global.radar.hydro.annual.anomaly=global.radar.hydro.annual[, , 1998-1993+1]-global.radar.hydro.annual.mean
amazon.radar.anomaly.1998=global.radar.hydro.annual.anomaly[401:520, 281:400]  # -20 
africa.radar.anomaly.1998=global.radar.hydro.annual.anomaly[741:860, 321:400]  # 10 190*4=760 start from 5E 

global.radar.hydro.annual.anomaly=global.radar.hydro.annual[, , 2016-1993+1]-global.radar.hydro.annual.mean
amazon.radar.anomaly.2016=global.radar.hydro.annual.anomaly[401:520, 281:400]  # -20 
africa.radar.anomaly.2016=global.radar.hydro.annual.anomaly[741:860, 321:400]  # start from 5E 185*4+1=741


load('CHIRPS_basins_mcwd_hydro_zscore_1982-2023_newbase1992-2020.Rdata')

amazon.mcwd.hydro.zscore[amazon.mcwd.hydro.zscore>0]=NA
amazon.mcwd.hydro.zscore[amazon.mcwd.hydro.zscore<(-3)]=-3


africa.mcwd.hydro.zscore[africa.mcwd.hydro.zscore>0]=NA
africa.mcwd.hydro.zscore[africa.mcwd.hydro.zscore<(-3)]=-3


amazon.mcwd.hydro.zscore.1998=amazon.mcwd.hydro.zscore[, , 1998-1982+1]
amazon.mcwd.hydro.zscore.2016=amazon.mcwd.hydro.zscore[, , 2016-1982+1]
amazon.mcwd.hydro.zscore.2023=amazon.mcwd.hydro.zscore[, , 2023-1982+1]

africa.mcwd.hydro.zscore.1998=africa.mcwd.hydro.zscore[, , 1998-1982+1]
africa.mcwd.hydro.zscore.2016=africa.mcwd.hydro.zscore[, , 2016-1982+1]
africa.mcwd.hydro.zscore.2023=africa.mcwd.hydro.zscore[, , 2023-1982+1]


amazon.mcwd.hydro.zscore.1998[amazon.tmf.005<0.95]=NA
amazon.mcwd.hydro.zscore.2016[amazon.tmf.005<0.95]=NA
amazon.mcwd.hydro.zscore.2023[amazon.tmf.005<0.95]=NA

africa.mcwd.hydro.zscore.1998[africa.tmf.005.match<0.95]=NA
africa.mcwd.hydro.zscore.2016[africa.tmf.005.match<0.95]=NA
africa.mcwd.hydro.zscore.2023[africa.tmf.005.match<0.95]=NA


load('/Users/yyao/Documents/futureAmazon/CHIRPS_hydro_conform_1982-2023_newbase1992-2020.Rdata')
#load(file='CHIRPS_hydro_conform_1982-2023_newbase1992-2020.Rdata')

amazon.hydro.conform.1998=amazon.hydro.conform[, , 1998-1982+1]
amazon.hydro.conform.2016=amazon.hydro.conform[, , 2016-1982+1]
amazon.hydro.conform.2023=amazon.hydro.conform[, , 2023-1982+1]

africa.hydro.conform.1998=africa.hydro.conform[, , 1998-1982+1]
africa.hydro.conform.2016=africa.hydro.conform[, , 2016-1982+1]
africa.hydro.conform.2023=africa.hydro.conform[, , 2023-1982+1]


amazon.hydro.conform.1998[amazon.tmf.005<0.95]=NA
amazon.hydro.conform.2016[amazon.tmf.005<0.95]=NA
amazon.hydro.conform.2023[amazon.tmf.005<0.95]=NA

amazon.hydro.conform.1998[is.na(amazon.mcwd.hydro.zscore.1998)]=NA
amazon.hydro.conform.2016[is.na(amazon.mcwd.hydro.zscore.2016)]=NA
amazon.hydro.conform.2023[is.na(amazon.mcwd.hydro.zscore.2023)]=NA


africa.hydro.conform.1998[africa.tmf.005.match<0.95]=NA
africa.hydro.conform.2016[africa.tmf.005.match<0.95]=NA
africa.hydro.conform.2023[africa.tmf.005.match<0.95]=NA

africa.hydro.conform.2023[is.na(africa.mcwd.hydro.zscore.2023)]=NA
africa.hydro.conform.2016[is.na(africa.mcwd.hydro.zscore.2016)]=NA
africa.hydro.conform.1998[is.na(africa.mcwd.hydro.zscore.1998)]=NA

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

amazon.radar.anomaly.1998.regrid=regrid_cscat_chirps_amazon(amazon.radar.anomaly.1998)
amazon.radar.anomaly.2016.regrid=regrid_cscat_chirps_amazon(amazon.radar.anomaly.2016)



projlat<-CRS("+proj=longlat +datum=WGS84 +no_defs")
amazon.Rmean<-raster(nrow=600,ncol=600)  # 0.05 deg 
extent(amazon.Rmean)<-c(-80, -50,-20,10)
amazon.Rmean[is.na(amazon.Rmean)]<-1
projection(amazon.Rmean)<-projlat
regrid_era5_chirps_amazon=function(mat) {
    #reproject1d<-function(mat,Rmean){
    pr_mask<-raster(t(mat[, 300:1]))
    extent(pr_mask)<-c(-80,-50,-20,10)
    projection(pr_mask)<-projlat
        
    meanpr_proj<-projectRaster(from = pr_mask,to=amazon.Rmean)
    mat_pr<-t(as.matrix(meanpr_proj))[,600:1]
    return(mat_pr)
    #}
}

amazon.tair.anomaly.1998=amazon.t2m.annual.anomaly[1:300, 1:300, 1998-1981+1]
amazon.tair.anomaly.2016=amazon.t2m.annual.anomaly[1:300, 1:300, 2016-1981+1]
amazon.tair.anomaly.1998.regrid=regrid_era5_chirps_amazon(amazon.tair.anomaly.1998)
amazon.tair.anomaly.2016.regrid=regrid_era5_chirps_amazon(amazon.tair.anomaly.2016)



congo.Rmean<-raster(nrow=400,ncol=600)  # 0.05 deg 
extent(congo.Rmean)<-c(5,35,-10,10)
congo.Rmean[is.na(congo.Rmean)]<-1
projection(congo.Rmean)<-projlat
regrid_era5_chirps_congo=function(mat) {
    #reproject1d<-function(mat,Rmean){
    pr_mask<-raster(t(mat[, 200:1]))
    extent(pr_mask)<-c(5,35,-10,10)
    projection(pr_mask)<-projlat
        
    meanpr_proj<-projectRaster(from = pr_mask,to=congo.Rmean)
    mat_pr<-t(as.matrix(meanpr_proj))[,400:1]
    return(mat_pr)
    #}
}
africa.tair.anomaly.1998=congo.t2m.annual.anomaly[51:350, 1:200, 1998-1981+1]
africa.tair.anomaly.2016=congo.t2m.annual.anomaly[51:350, 1:200, 2016-1981+1]
congo.tair.anomaly.1998.regrid=regrid_era5_chirps_congo(africa.tair.anomaly.1998)
congo.tair.anomaly.2016.regrid=regrid_era5_chirps_congo(africa.tair.anomaly.2016)



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


congo.radar.anomaly.1998.regrid=regrid_cscat_chirps_congo(africa.radar.anomaly.1998)
congo.radar.anomaly.2016.regrid=regrid_cscat_chirps_congo(africa.radar.anomaly.2016)

###############



count_area=function(radar.anomaly, hydro.confirm, mcwd.zscore, forest.mask,forest.area)  {
    radar.anomaly=radar.anomaly*(forest.mask>=0.95)
    radar.anomaly[radar.anomaly==0]=NA
    mcwd.zscore[mcwd.zscore<(-3)]=-3
    bins.1=c(-1, -1.645, -1.96, -2.576)
    bins.2=c(-1.645, -1.96, -2.576, -3)
    bk.bins.area=rep(NA, 16)
    dim(bk.bins.area)=c(2,4, 2)
    for (bb in 1:4) {
        sbk=(mcwd.zscore>=bins.2[bb]) * (mcwd.zscore<bins.1[bb]) * (hydro.confirm>0) * (forest.mask>=0.95) *(radar.anomaly)*forest.area
        sbk[sbk==0]=NA
        sarea=(!is.na(sbk))*forest.area
        bk.bins.area[1, bb, 1]=sum(sbk, na.rm=TRUE)/sum(sarea, na.rm=TRUE)  # area-weighted values 
        bk.bins.area[1, bb, 2]=sum(sarea, na.rm=TRUE)/sum((forest.mask>=0.95)*forest.area, na.rm=TRUE)*100
    }
    for (bb in 1:4) {
        sbk=(mcwd.zscore>=bins.2[bb]) * (mcwd.zscore<bins.1[bb]) * (hydro.confirm<0) * (forest.mask>=0.95) *(radar.anomaly)*forest.area
        sbk[sbk==0]=NA
        sarea=(!is.na(sbk))*forest.area
        bk.bins.area[2, bb, 1]=sum(sbk, na.rm=TRUE)/sum(sarea, na.rm=TRUE)  # area-weighted values 
        bk.bins.area[2, bb, 2]=sum(sarea, na.rm=TRUE)/sum((forest.mask>=0.95)*forest.area, na.rm=TRUE)*100
    }
    return(bk.bins.area)
}

amazon.1998.count.area=count_area(amazon.radar.anomaly.1998.regrid, amazon.hydro.conform.1998, amazon.mcwd.hydro.zscore.1998, amazon.tmf.005, amazon.area) 
amazon.2016.count.area=count_area(amazon.radar.anomaly.2016.regrid, amazon.hydro.conform.2016, amazon.mcwd.hydro.zscore.2016, amazon.tmf.005, amazon.area) 

congo.1998.count.area=count_area(congo.radar.anomaly.1998.regrid, africa.hydro.conform.1998, africa.mcwd.hydro.zscore.1998, africa.tmf.005.match, congo.area) 
congo.2016.count.area=count_area(congo.radar.anomaly.2016.regrid, africa.hydro.conform.2016, africa.mcwd.hydro.zscore.2016, africa.tmf.005.match, congo.area) 


plot_count_area_wet=function(count_bk_area, title) {
    plot(1:4, 1:4, col='white', xlim=c(-50,50), ylim=c(-0.15, 0.04),main=title, cex.main=2, las=1,
    xlab='Area fraction (%)', ylab='', cex.axis=2, cex.lab=2, xaxt='n')

    mtext(side=2, text=expression(Delta~"BK (dB)"), line=5, cex.axis=2, cex=2)

    bk.value=count_bk_area[1, , 1]
    bk.area=count_bk_area[1, , 2]  # wet 
    ll=0
    for (bb in 1:4) {
        rr=ll+bk.area[bb]
        rect(ll, 0, rr, bk.value[bb], col=colos[bb])
        ll=ll+bk.area[bb]
    }

    bk.value=count_bk_area[2, , 1]  # dry 
    bk.area=count_bk_area[2, , 2]
    ll=0
    for (bb in 1:4) {
        rr=ll-bk.area[bb]
        rect(ll, 0, rr, bk.value[bb], col=colos[bb])
        ll=ll-bk.area[bb]
    }

    text(-29, 0.02, 'Dry season drought', cex=2)
    text(29, 0.02, 'Wet season drought', cex=2)
}

tiff('figure3_Oct17.tiff', width=1100, height=600)
#colos=brewer.pal(5,'Reds')[2:5]
colos=c('#fed976','#fd8d3c','#e31a1c','#800026')
par(mfrow=c(2,2))
par(mar=c(5,7.5,3,3))
plot_count_area_wet(amazon.1998.count.area, '(a) Amazon 1997/98')
abline(h=0, lty=2, lwd=2)
abline(v=0, lty=2, lwd=2)
legend('bottomright', legend=c('Mild','Moderate','Severe','Extreme'), text.col=colos, bty='n', cex=2)
axis(1, at=seq(-40,40,20), labels=c('40','20','0','20','40'), cex=2, cex.axis=2)

par(mar=c(5,7.5,3,3))
plot_count_area_wet(amazon.2016.count.area, '(b) Amazon 2015/16')
abline(h=0, lty=2, lwd=2)
abline(v=0, lty=2, lwd=2)
axis(1, at=seq(-40,40,20), labels=c('40','20','0','20','40'), cex=2, cex.axis=2)

par(mar=c(5,7.5,3,3))
plot_count_area_wet(congo.1998.count.area, '(c) C. Africa 1997/98')
abline(h=0, lty=2, lwd=2)
abline(v=0, lty=2, lwd=2)
axis(1, at=seq(-40,40,20), labels=c('40','20','0','20','40'), cex=2, cex.axis=2)

par(mar=c(5,7.5,3,3))
plot_count_area_wet(congo.2016.count.area, '(d) C. Africa 2015/16')
abline(h=0, lty=2, lwd=2)
abline(v=0, lty=2, lwd=2)
axis(1, at=seq(-40,40,20), labels=c('40','20','0','20','40'), cex=2, cex.axis=2)

dev.off()
