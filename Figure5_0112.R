
load(file='/Users/yyao/Downloads/Cscat/0127radar_global_grid_hydro_annual_1992-2022.Rdata')
global.radar.hydro.annual.mean=apply(global.radar.hydro.annual[, , (1993-1993+1):(2020-1993+1)], c(1,2), mean, na.rm=TRUE) 

global.radar.hydro.annual.anomaly=global.radar.hydro.annual[, , 1998-1993+1]-global.radar.hydro.annual.mean
amazon.radar.anomaly.1998=global.radar.hydro.annual.anomaly[401:520, 281:400]  # -20 
africa.radar.anomaly.1998=global.radar.hydro.annual.anomaly[741:860, 321:400]  # 10 190*4=760 start from 5E 

global.radar.hydro.annual.anomaly=global.radar.hydro.annual[, , 2016-1993+1]-global.radar.hydro.annual.mean
amazon.radar.anomaly.2016=global.radar.hydro.annual.anomaly[401:520, 281:400]  # -20 
africa.radar.anomaly.2016=global.radar.hydro.annual.anomaly[741:860, 321:400]  # start from 5E 185*4+1=741



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


amazon.radar.anomaly.2016.regrid=regrid_cscat_chirps_amazon(amazon.radar.anomaly.2016)

congo.radar.anomaly.2016.regrid=regrid_cscat_chirps_congo(africa.radar.anomaly.2016)



load(file='amazon_t2m_hydro_anomaly_oct.Rdata')
load(file='congo_t2m_hydro_anomaly_oct.Rdata')

amazon.tair.anomaly.1998=letter.t2m.annual.oct.anomaly.1998[, 300:1]
amazon.tair.anomaly.2016=letter.t2m.annual.oct.anomaly.2016[, 300:1]
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

africa.tair.anomaly.1998=congo.t2m.annual.oct.anomaly.1998[51:350, 200:1]
africa.tair.anomaly.2016=congo.t2m.annual.oct.anomaly.2016[51:350, 200:1]
congo.tair.anomaly.1998.regrid=regrid_era5_chirps_congo(africa.tair.anomaly.1998)
congo.tair.anomaly.2016.regrid=regrid_era5_chirps_congo(africa.tair.anomaly.2016)

load(file='CHIRPS_basins_mcwd_hydro_anomaly_1992-2020_newbase1992-2020.Rdata')
# amazon.mcwd.hydro.tao.anomaly, africa.mcwd.hydro.tao.anomaly
####
#amazon.1998.mcwd.ctl=amazon.annual.mcwd.sens.2[,, 1]*amazon.mcwd.annual.anomaly.tao[,,1998-1992+1]/100
#amazon.1998.tair.ctl=amazon.annual.tair.sens.2[,, 1]*amazon.t2m.annual.anomaly.tao.deg[,,1998-1992+1]

amazon.2016.mcwd.ctl=amazon.annual.mcwd.sens.2[,, 2]*amazon.mcwd.hydro.tao.anomaly[, , 2016-1992+1]/100 #amazon.mcwd.annual.anomaly.tao[,,2016-1992+1]/100
amazon.2016.tair.ctl=amazon.annual.tair.sens.2[,, 2]*amazon.tair.anomaly.2016.regrid  #amazon.t2m.annual.anomaly.tao.deg[,,2016-1992+1]

amazon.2016.mcwd.sens=amazon.annual.mcwd.sens.2[,, 1]*amazon.mcwd.hydro.tao.anomaly[, , 2016-1992+1]/100  #amazon.mcwd.annual.anomaly.tao[,,2016-1992+1]/100
amazon.2016.tair.sens=amazon.annual.tair.sens.2[,, 1]*amazon.tair.anomaly.2016.regrid  #amazon.t2m.annual.anomaly.tao.deg[,,2016-1992+1]

amazon.2016.mcwd.force=amazon.annual.mcwd.sens.2[,, 2]*amazon.mcwd.hydro.tao.anomaly[, , 1998-1992+1]/100  #amazon.mcwd.annual.anomaly.tao[,,1998-1992+1]/100
amazon.2016.tair.force=amazon.annual.tair.sens.2[,, 2]*amazon.tair.anomaly.1998.regrid  #amazon.t2m.annual.anomaly.tao.deg[,,1998-1992+1]

#################
y0=amazon.radar.anomaly.2016.regrid  #amazon.radar.tao.anomaly.deg.annual[,,2016-1992+1]
y1=amazon.2016.mcwd.ctl+amazon.2016.tair.ctl
y2=amazon.2016.mcwd.sens+amazon.2016.tair.ctl
y3=amazon.2016.mcwd.force+amazon.2016.tair.ctl

y4=amazon.2016.mcwd.ctl+amazon.2016.tair.sens
y5=amazon.2016.mcwd.ctl+amazon.2016.tair.force

amazon.y=array(NA, dim=c(600,600,6))
amazon.y[, , 1]=y0
amazon.y[,,2]=y1
amazon.y[,,3]=y2
amazon.y[,,4]=y3
amazon.y[,,5]=y4
amazon.y[,,6]=y5

amazon.region.test=array(NA, dim=c(6,4))
for (rr in 1:4) {
    for (tt in 1:6) {
        sx=amazon.y[,,tt]*(amazon.tmf.005>=0.95)*(amazon.mask5==rr)
        sx[sx==0]=NA
        amazon.region.test[tt, rr]=mean(sx, na.rm=TRUE)
    }
}



#congo.1998.mcwd.ctl=congo.annual.mcwd.sens.2[,, 1]*africa.mcwd.annual.anomaly.tao[,,1998-1992+1]/100
#congo.1998.tair.ctl=congo.annual.tair.sens.2[,, 1]*congo.t2m.annual.anomaly.tao.deg[,,1998-1992+1]

congo.2016.mcwd.ctl=congo.annual.mcwd.sens.2[,, 2]*africa.mcwd.hydro.tao.anomaly[, , 2016-1992+1]/100   #africa.mcwd.annual.anomaly.tao[,,2016-1992+1]/100
congo.2016.tair.ctl=congo.annual.tair.sens.2[,, 2]*congo.tair.anomaly.2016.regrid   #congo.t2m.annual.anomaly.tao.deg[,,2016-1992+1]

congo.2016.mcwd.sens=congo.annual.mcwd.sens.2[,, 1]*africa.mcwd.hydro.tao.anomaly[, , 2016-1992+1]/100    #africa.mcwd.annual.anomaly.tao[,,2016-1992+1]/100
congo.2016.tair.sens=congo.annual.tair.sens.2[,, 1]*congo.tair.anomaly.2016.regrid  #congo.t2m.annual.anomaly.tao.deg[,,2016-1992+1]

congo.2016.mcwd.force=congo.annual.mcwd.sens.2[,, 2]*africa.mcwd.hydro.tao.anomaly[, , 1998-1992+1]/100    #africa.mcwd.annual.anomaly.tao[,,1998-1992+1]/100
congo.2016.tair.force=congo.annual.tair.sens.2[,, 2]*congo.tair.anomaly.1998.regrid  #congo.t2m.annual.anomaly.tao.deg[,,1998-1992+1]

#################
y0=congo.radar.anomaly.2016.regrid  #africa.radar.tao.anomaly.deg.annual[, , 2016-1992+1]
y1=congo.2016.mcwd.ctl+congo.2016.tair.ctl
y2=congo.2016.mcwd.sens+congo.2016.tair.ctl
y3=congo.2016.mcwd.force+congo.2016.tair.ctl

y4=congo.2016.mcwd.ctl+congo.2016.tair.sens
y5=congo.2016.mcwd.ctl+congo.2016.tair.force

congo.y=array(NA, dim=c(600,400,6))
congo.y[, , 1]=y0
congo.y[,,2]=y1
congo.y[,,3]=y2
congo.y[,,4]=y3
congo.y[,,5]=y4
congo.y[, , 6]=y5

congo.region.test=array(NA, dim=c(6,4))
for (rr in 1:4) {
    for (tt in 1:6) {
        sx=congo.y[,,tt]*(africa.tmf.005.match>=0.95)*(africa.mask5==rr)
        sx[sx==0]=NA
        congo.region.test[tt, rr]=mean(sx, na.rm=TRUE)
    }
}

barplot(congo.region.test, beside=TRUE)


tiff('figure5_factor_0112.tiff', width=800, height=700)
colos=c('#d9d9d9','#fccde5','#80b1d3','#8dd3c7','#fb8072','#fdb462')
par(mfrow=c(2,1))

par(mar=c(5,7,3,3))
barplot(amazon.region.test, beside=TRUE, col=colos, las=1, 
names.arg=c(expression('R'[mild]^'dry'),expression('R'[mild]^'wet'),expression('R'[severe]^'dry'),expression('R'[severe]^'wet')), 
cex.names=2,  
main='(a) Amazon', cex.main=2, ylab='',cex.lab=2, cex.axis=2)
mtext(side=2, expression(Delta~'BK (dB)'), line=5, cex=2)

par(mar=c(5,7,3,3))
barplot(congo.region.test, beside=TRUE, col=colos, las=1, 
names.arg=c(expression('R'[mild]^'dry'),expression('R'[mild]^'wet'),expression('R'[severe]^'dry'),expression('R'[severe]^'wet')), 
cex.names=2, 
main='(b) C. Africa', cex.main=2, ylab='',cex.lab=2, cex.axis=2)
mtext(side=2, expression(Delta~'BK (dB)'), line=5, cex=2)

t1=expression('T1: held'~gamma[MCWD])
t2=expression('T2: held'~Delta~'MCWD')
t3=expression('T3: held'~gamma[Tair])
t4=expression('T4: held'~Delta~'Tair')
legend('bottomleft', legend=c('Obs','T0: MLR',
t1, t2, t3, t4), text.col=colos, cex=2, bty='n')

dev.off()

