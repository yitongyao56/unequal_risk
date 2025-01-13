load(file='/users/yyao/Documents/futureAmazon/amazon_fraction_tmf_deg01.Rdata')
load(file='/users/yyao/Documents/futureAmazon/amazon_fraction_tmf_deg005.Rdata')

#load(file='/users/yyao/Documents/futureAmazon/africa_fraction_tmf_deg05.Rdata')
load(file='/users/yyao/Documents/futureAmazon/africa_fraction_tmf_deg01_June22.Rdata')
load(file='/users/yyao/Documents/futureAmazon/africa_fraction_tmf_deg005_June22.Rdata')

# 10 to 35 -10 to 10 
#africa.tmf.05.match=africa.tmf.05[1:50, ]
africa.tmf.01.match=africa.tmf.01[51:350, ]
africa.tmf.005.match=africa.tmf.005[101:700, ]  # start from 5 

library(raster)
library(fields)
library(maps)

figure5_barplot=function(sensitivity, forest.tmf, s1, s2, s3, s4, forest.area) {
    out=rep(NA, 4)
    sensitivity[forest.tmf<0.95]=NA
    
    draw.trend=sensitivity*s1
    draw.trend[draw.trend==0]=NA
    draw.trend[draw.trend>1]=1
    draw.trend[draw.trend<(-1)]=-1
    strend=draw.trend*forest.area
    sarea=(!is.na(draw.trend))*forest.area
    out[1]=sum(strend, na.rm=TRUE)/sum(sarea, na.rm=TRUE)

    draw.trend=sensitivity*s2
    draw.trend[draw.trend==0]=NA
    draw.trend[draw.trend>1]=1
    draw.trend[draw.trend<(-1)]=-1
    strend=draw.trend*forest.area
    sarea=(!is.na(draw.trend))*forest.area
    out[2]=sum(strend, na.rm=TRUE)/sum(sarea, na.rm=TRUE)

    draw.trend=sensitivity*s3
    draw.trend[draw.trend==0]=NA
    draw.trend[draw.trend>1]=1
    draw.trend[draw.trend<(-1)]=-1
    strend=draw.trend*forest.area
    sarea=(!is.na(draw.trend))*forest.area
    out[3]=sum(strend, na.rm=TRUE)/sum(sarea, na.rm=TRUE)

    draw.trend=sensitivity*s4
    draw.trend[draw.trend==0]=NA
    draw.trend[draw.trend>1]=1
    draw.trend[draw.trend<(-1)]=-1
    strend=draw.trend*forest.area
    sarea=(!is.na(draw.trend))*forest.area
    out[4]=sum(strend, na.rm=TRUE)/sum(sarea, na.rm=TRUE)

    return(out)
}



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




load(file='/Users/yyao/Documents/futureAmazon/ERA5_amazon_ssrd_annual.Rdata')
Amazon.basin.ssrd.annual.base=apply(Amazon.basin.ssrd.annual[, , 12:40], c(1,2), mean, na.rm=TRUE) # 2000-1981+1=20 1992-1981+1=12 2020-1981+1=40 
Amazon.basin.ssrd.annual.base[Amazon.basin.ssrd.annual.base==0]=NA
Amazon.basin.ssrd.annual.anomaly=array(NA, dim=c(300,300,43))
for (year in 1:43) {
    Amazon.basin.ssrd.annual.anomaly[, , year]=Amazon.basin.ssrd.annual[,,year]-Amazon.basin.ssrd.annual.base
}

load(file='/Users/yyao/Documents/futureAmazon/ERA5_amazon_vpd_annual.Rdata')
amazon.vpd.annual.base=apply(amazon.vpd.annual[, , 12:40], c(1,2), mean, na.rm=TRUE) # 2000-1981+1=20 1992-1981+1=12 2020-1981+1=40 
amazon.vpd.annual.base[amazon.vpd.annual.base==0]=NA
amazon.vpd.annual.anomaly=array(NA, dim=c(301,301,43))
for (year in 1:43) {
    amazon.vpd.annual.anomaly[, , year]=amazon.vpd.annual[,,year]-amazon.vpd.annual.base
}

load(file='/Users/yyao/Documents/futureAmazon/ERA5_amazon_t2m_annual.Rdata')
amazon.t2m.annual.base=apply(amazon.t2m.annual[, , 12:40], c(1,2), mean, na.rm=TRUE)
amazon.t2m.annual.base[amazon.t2m.annual.base==0]=NA
amazon.t2m.annual.anomaly=array(NA, dim=c(301,301,43))
for (year in 1:43) {
    amazon.t2m.annual.anomaly[, , year]=amazon.t2m.annual[,,year]-amazon.t2m.annual.base
}

amazon.vpd.annual.anomaly.tao=amazon.vpd.annual.anomaly[1:300, 1:300, (1992-1981+1):(2020-1981+1)]
amazon.t2m.annual.anomaly.tao=amazon.t2m.annual.anomaly[1:300, 1:300, (1992-1981+1):(2020-1981+1)]
Amazon.basin.ssrd.annual.anomaly.tao=Amazon.basin.ssrd.annual.anomaly[1:300, 1:300, (1992-1981+1):(2020-1981+1)]

amazon.vpd.annual.anomaly.tao.deg=array(NA, dim=c(600,600,29))
amazon.t2m.annual.anomaly.tao.deg=array(NA, dim=c(600,600,29))
amazon.ssrd.annual.anomaly.tao.deg=array(NA, dim=c(600,600,29))
for (year in 1:29) {
    amazon.vpd.annual.anomaly.tao.deg[, , year]=regrid_era5_chirps_amazon(amazon.vpd.annual.anomaly.tao[, , year])
    amazon.t2m.annual.anomaly.tao.deg[, , year]=regrid_era5_chirps_amazon(amazon.t2m.annual.anomaly.tao[, , year])
    amazon.ssrd.annual.anomaly.tao.deg[, , year]=regrid_era5_chirps_amazon(Amazon.basin.ssrd.annual.anomaly.tao[, , year])
}



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


load(file='/Users/yyao/Documents/futureAmazon/ERA5_africa_vpd_annual.Rdata')
congo.vpd.annual.base=apply(Africa.vpd.annual[, , 12:40], c(1,2), mean, na.rm=TRUE)
congo.vpd.annual.base[congo.vpd.annual.base==0]=NA
congo.vpd.annual.anomaly=array(NA, dim=c(401,201,43))
for (year in 1:43) {
    congo.vpd.annual.anomaly[, , year]=Africa.vpd.annual[,,year]-congo.vpd.annual.base
}

load(file='/Users/yyao/Documents/futureAmazon/ERA5_africa_t2m_annual.Rdata')
congo.t2m.annual.base=apply(Africa.t2m.annual[, , 12:40], c(1,2), mean, na.rm=TRUE)
congo.t2m.annual.base[congo.t2m.annual.base==0]=NA
congo.t2m.annual.anomaly=array(NA, dim=c(401,201,43))
for (year in 1:43) {
    congo.t2m.annual.anomaly[, , year]=Africa.t2m.annual[,,year]-congo.t2m.annual.base
}


load(file='/Users/yyao/Documents/futureAmazon/ERA5_africa_ssrd_annual.Rdata')
congo.ssrd.annual.base=apply(Africa.ssrd.annual[, , 12:40], c(1,2), mean, na.rm=TRUE)
congo.ssrd.annual.base[congo.ssrd.annual.base==0]=NA
congo.ssrd.annual.anomaly=array(NA, dim=c(401,201,43))
for (year in 1:43) {
    congo.ssrd.annual.anomaly[, , year]=Africa.ssrd.annual[,,year]-congo.ssrd.annual.base
}


congo.vpd.annual.anomaly.tao=congo.vpd.annual.anomaly[51:350, 1:200, (1992-1981+1):(2020-1981+1)]
congo.t2m.annual.anomaly.tao=congo.t2m.annual.anomaly[51:350, 1:200, (1992-1981+1):(2020-1981+1)]
congo.ssrd.annual.anomaly.tao=congo.ssrd.annual.anomaly[51:350, 1:200, (1992-1981+1):(2020-1981+1)]

congo.vpd.annual.anomaly.tao.deg=array(NA, dim=c(600,400,29))
congo.t2m.annual.anomaly.tao.deg=array(NA, dim=c(600,400,29))
congo.ssrd.annual.anomaly.tao.deg=array(NA, dim=c(600,400,29))
for (year in 1:29) {
    congo.vpd.annual.anomaly.tao.deg[, , year]=regrid_era5_chirps_congo(congo.vpd.annual.anomaly.tao[, , year])
    congo.t2m.annual.anomaly.tao.deg[, , year]=regrid_era5_chirps_congo(congo.t2m.annual.anomaly.tao[, , year])
    congo.ssrd.annual.anomaly.tao.deg[, , year]=regrid_era5_chirps_congo(congo.ssrd.annual.anomaly.tao[, , year])
}




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

load(file='/Users/yyao/Documents/futureAmazon/Amazon_radar_tao_annual_anomaly_025deg_1992-2020.Rdata')
amazon.radar.tao.anomaly.deg.annual=array(NA, dim=c(600,600,29))
for (year in 1:29) {
    amazon.radar.tao.anomaly.deg.annual[, , year]=regrid_cscat_chirps_amazon(amazon.radar.tao.annual.anomaly[, , year])
}
amazon.radar.tao.anomaly.deg.annual[amazon.radar.tao.anomaly.deg.annual==0]=NA


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
load(file='Africa_radar_tao_annual_anomaly_025deg_1992-2020.Rdata')
africa.radar.tao.anomaly.deg.annual=array(NA, dim=c(600,400,29))
for (year in 1:29) {
    africa.radar.tao.anomaly.deg.annual[, , year]=regrid_cscat_chirps_congo(africa.radar.tao.annual.anomaly[, , year])
}
africa.radar.tao.anomaly.deg.annual[africa.radar.tao.anomaly.deg.annual==0]=NA


load(file='chirps_mcwd_annual_anomaly_1981=2023_newbase1992-2020.Rdata')
amazon.mcwd.annual.anomaly.tao=amazon.mcwd.annual.anomaly[, , (1992-1981+1):(2020-1981+1)]
africa.mcwd.annual.anomaly.tao=africa.mcwd.annual.anomaly[, , (1992-1981+1):(2020-1981+1)]


amazon.annual.mcwd.sens.2=array(NA, dim=c(600,600,2))
amazon.annual.tair.sens.2=array(NA, dim=c(600,600,2))
amazon.annual.deg.r2.2=array(NA, dim=c(600,600,2))

amazon.annual.anoma.sens.2=array(NA, dim=c(600,600,2))
amazon.annual.anoma.tair.2=array(NA, dim=c(600,600,2))

amazon.annual.anoma.cor.2=array(NA, dim=c(600,600,2))

bb=1
for (rr in 1:600 )  {
    for (cc in 1:600) {
        y1=amazon.radar.tao.anomaly.deg.annual[rr, cc, c(1:15)]
        x1=amazon.mcwd.annual.anomaly.tao[rr, cc, c(1:15)]
        x2=amazon.vpd.annual.anomaly.tao.deg[rr, cc, c(1:15)]
        x3=amazon.t2m.annual.anomaly.tao.deg[rr, cc, c(1:15)]
        if (sum(!is.na(y1))>=8 & sum(!is.na(x1))>=8 & sum(!is.na(x2))>=8 & sum(!is.na(x3))>=8) {
            lmm=lm(y1~x1+x3)
            amazon.annual.mcwd.sens.2[rr, cc, bb]=lmm$coefficients[2]*100  # change of radar signal anomaly per 100mm increase in MCWD 
            amazon.annual.tair.sens.2[rr, cc, bb]=lmm$coefficients[3]
            xx=summary(lmm)
            amazon.annual.deg.r2.2[rr, cc, bb]=xx$r.squared

            amazon.annual.anoma.sens.2[rr, cc, bb]=amazon.annual.mcwd.sens.2[rr, cc, bb]*amazon.mcwd.annual.anomaly.tao[rr, cc,1998-1992+1]/100
            #amazon.annual.anoma.tair.2[rr, cc, bb]=amazon.annual.tair.sens.2[rr, cc, bb]*amazon.vpd.annual.anomaly.tao.deg[rr, cc,1998-1992+1]
            amazon.annual.anoma.tair.2[rr, cc, bb]=amazon.annual.tair.sens.2[rr, cc, bb]*amazon.t2m.annual.anomaly.tao.deg[rr, cc,1998-1992+1]

            newy=x1*lmm$coefficients[2]+x3*lmm$coefficients[3]
            amazon.annual.anoma.cor.2[rr, cc, bb]=cor(y1, newy)
        }
    }
}
bb=2
for (rr in 1:600 )  {
    for (cc in 1:600) {
        y1=amazon.radar.tao.anomaly.deg.annual[rr, cc, c(15:29)]
        x1=amazon.mcwd.annual.anomaly.tao[rr, cc, c(15:29)]
        x2=amazon.vpd.annual.anomaly.tao.deg[rr, cc, c(15:29)]
        x3=amazon.t2m.annual.anomaly.tao.deg[rr, cc, c(15:29)]
        if (sum(!is.na(y1))>=8 & sum(!is.na(x1))>=8 & sum(!is.na(x2))>=8 & sum(!is.na(x3))>=8) {
            lmm=lm(y1~x1+x3)
            amazon.annual.mcwd.sens.2[rr, cc, bb]=lmm$coefficients[2]*100  # change of radar signal anomaly per 100mm increase in MCWD 
            amazon.annual.tair.sens.2[rr, cc, bb]=lmm$coefficients[3]
            #amazon.annual.tair.sens.2[rr, cc, bb]=lmm$coefficients[4]
            xx=summary(lmm)
            amazon.annual.deg.r2.2[rr, cc, bb]=xx$r.squared

            amazon.annual.anoma.sens.2[rr, cc, bb]=amazon.annual.mcwd.sens.2[rr, cc, bb]*amazon.mcwd.annual.anomaly.tao[rr, cc,2016-1992+1]/100
            #amazon.annual.anoma.vpd.2[rr, cc, bb]=amazon.annual.deg.vpd.2[rr, cc, bb]*amazon.vpd.annual.anomaly.tao.deg[rr, cc,2016-1992+1]
            amazon.annual.anoma.tair.2[rr, cc, bb]=amazon.annual.tair.sens.2[rr, cc, bb]*amazon.t2m.annual.anomaly.tao.deg[rr, cc,2016-1992+1]

            newy=x1*lmm$coefficients[2]+x3*lmm$coefficients[3]
            amazon.annual.anoma.cor.2[rr, cc, bb]=cor(y1, newy)
        }
    }
}



congo.annual.mcwd.sens.2=array(NA, dim=c(600,400,2))
congo.annual.tair.sens.2=array(NA, dim=c(600,400,2))
congo.annual.deg.vpd.2=array(NA, dim=c(600,400,2))
congo.annual.deg.cor2=array(NA, dim=c(600,400,2))

congo.annual.anoma.sens.2=array(NA, dim=c(600,400,2))
congo.annual.anoma.tair.2=array(NA, dim=c(600,400,2))
congo.annual.anoma.vpd.2=array(NA, dim=c(600,400,2))

congo.annual.anoma.cor.2=array(NA, dim=c(600,400,2))

bb=1
for (rr in 1:600 )  {
    for (cc in 1:400) {
        y1=africa.radar.tao.anomaly.deg.annual[rr, cc, c(1:15)]
        x1=africa.mcwd.annual.anomaly.tao[rr, cc, c(1:15)]
        x2=congo.vpd.annual.anomaly.tao.deg[rr, cc, c(1:15)]
        x3=congo.t2m.annual.anomaly.tao.deg[rr, cc, c(1:15)]
        if (sum(!is.na(y1))>=8 & sum(!is.na(x1))>=8 & sum(!is.na(x2))>=8 & sum(!is.na(x3))>=8) {
            lmm=lm(y1~x1+x3)
            congo.annual.mcwd.sens.2[rr, cc, bb]=lmm$coefficients[2]*100  # change of radar signal anomaly per 100mm increase in MCWD 
            congo.annual.tair.sens.2[rr, cc, bb]=lmm$coefficients[3]
            #congo.annual.tair.sens.2[rr, cc, bb]=lmm$coefficients[4]
            xx=summary(lmm)
            #congo.annual.deg.r2.2[rr, cc, bb]=xx$r.squared

            congo.annual.anoma.sens.2[rr, cc, bb]=congo.annual.mcwd.sens.2[rr, cc, bb]*africa.mcwd.annual.anomaly.tao[rr, cc,1998-1992+1]/100
            #congo.annual.anoma.vpd.2[rr, cc, bb]=congo.annual.deg.vpd.2[rr, cc, bb]*congo.vpd.annual.anomaly.tao.deg[rr, cc,1998-1992+1]
            congo.annual.anoma.tair.2[rr, cc, bb]=congo.annual.tair.sens.2[rr, cc, bb]*congo.t2m.annual.anomaly.tao.deg[rr, cc,1998-1992+1]

            newy=x1*lmm$coefficients[2]+x3*lmm$coefficients[3]
            congo.annual.anoma.cor.2[rr, cc, bb]=cor(y1, newy)
        }
    }
}
bb=2
for (rr in 1:600 )  {
    for (cc in 1:400) {
        y1=africa.radar.tao.anomaly.deg.annual[rr, cc, c(15:29)]
        x1=africa.mcwd.annual.anomaly.tao[rr, cc, c(15:29)]
        x2=congo.vpd.annual.anomaly.tao.deg[rr, cc, c(15:29)]
        x3=congo.t2m.annual.anomaly.tao.deg[rr, cc, c(15:29)]
        if (sum(!is.na(y1))>=8 & sum(!is.na(x1))>=8 & sum(!is.na(x2))>=8 & sum(!is.na(x3))>=8) {
            lmm=lm(y1~x1+x3)
            congo.annual.mcwd.sens.2[rr, cc, bb]=lmm$coefficients[2]*100  # change of radar signal anomaly per 100mm increase in MCWD 
            congo.annual.tair.sens.2[rr, cc, bb]=lmm$coefficients[3]
            #congo.annual.tair.sens.2[rr, cc, bb]=lmm$coefficients[4]
            xx=summary(lmm)
            #congo.annual.deg.r2.2[rr, cc, bb]=xx$r.squared

            congo.annual.anoma.sens.2[rr, cc, bb]=congo.annual.mcwd.sens.2[rr, cc, bb]*africa.mcwd.annual.anomaly.tao[rr, cc,2016-1992+1]/100
            #congo.annual.anoma.vpd.2[rr, cc, bb]=congo.annual.deg.vpd.2[rr, cc, bb]*congo.vpd.annual.anomaly.tao.deg[rr, cc,1998-1992+1]
            congo.annual.anoma.tair.2[rr, cc, bb]=congo.annual.tair.sens.2[rr, cc, bb]*congo.t2m.annual.anomaly.tao.deg[rr, cc,2016-1992+1]

            #newy=x1*lmm$coefficients[2]+x2*lmm$coefficients[3]+x3*lmm$coefficients[4]
            #congo.annual.deg.cor2[rr, cc, bb]=cor(y1, newy)
            newy=x1*lmm$coefficients[2]+x3*lmm$coefficients[3]
            congo.annual.anoma.cor.2[rr, cc, bb]=cor(y1, newy)
        }
    }
}




load('CHIRPS_basins_mcwd_hydro_zscore_1982-2023_newbase1992-2020.Rdata')

#amazon.mcwd.hydro.zscore[amazon.mcwd.hydro.zscore>0]=NA
amazon.mcwd.hydro.zscore[amazon.mcwd.hydro.zscore<(-3)]=-3


#africa.mcwd.hydro.zscore[africa.mcwd.hydro.zscore>0]=NA
africa.mcwd.hydro.zscore[africa.mcwd.hydro.zscore<(-3)]=-3


amazon.mcwd.hydro.zscore.1998=amazon.mcwd.hydro.zscore[, , 1998-1982+1]  # 1982-2023 
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

#amazon.hydro.conform.1998[is.na(amazon.mcwd.hydro.zscore.1998)]=NA
#amazon.hydro.conform.2016[is.na(amazon.mcwd.hydro.zscore.2016)]=NA
#amazon.hydro.conform.2023[is.na(amazon.mcwd.hydro.zscore.2023)]=NA


africa.hydro.conform.1998[africa.tmf.005.match<0.95]=NA
africa.hydro.conform.2016[africa.tmf.005.match<0.95]=NA
africa.hydro.conform.2023[africa.tmf.005.match<0.95]=NA
#africa.hydro.conform.2023[is.na(africa.mcwd.hydro.zscore.2023)]=NA
#africa.hydro.conform.2016[is.na(africa.mcwd.hydro.zscore.2016)]=NA
#africa.hydro.conform.1998[is.na(africa.mcwd.hydro.zscore.1998)]=NA


#amazon_chirps_annual_mcwd_june_zscore_2016[amazon_chirps_annual_mcwd_june_zscore_2016>=0]=NA
amazon_chirps_select_zscore_2016=amazon.mcwd.hydro.zscore.2016 #amazon_chirps_annual_mcwd_june_zscore_2016
#amazon_chirps_select_zscore_2016[amazon_chirps_select_zscore_2016>=0]=NA
#amazon_chirps_select_zscore_2016[amazon_chirps_select_zscore_2016>0]=NA

amazon.s1=amazon_chirps_select_zscore_2016>=(-1) & (amazon.hydro.conform.2016<0)  # non severe dry 
amazon.s2=amazon_chirps_select_zscore_2016>=(-1) & (amazon.hydro.conform.2016>0)  # non severe wet 
amazon.s3=amazon_chirps_select_zscore_2016<(-1) & (amazon.hydro.conform.2016<0)  # severe dry 
amazon.s4=amazon_chirps_select_zscore_2016<(-1) & (amazon.hydro.conform.2016>0)  # severe wet 

amazon.mask5=array(NA, dim=c(600,600))
amazon.mask5[amazon.s1==1]=1
amazon.mask5[amazon.s2==1]=2
amazon.mask5[amazon.s3==1]=3
amazon.mask5[amazon.s4==1]=4
amazon.mask5=amazon.mask5*(amazon.tmf.005>=0.95)
amazon.mask5[amazon.mask5==0]=NA

load('area_amazon_congo_005.Rdata')

amazon.mcwd.bar.1=figure5_barplot(amazon.annual.mcwd.sens.2[,,1], amazon.tmf.005, amazon.s1, amazon.s2, amazon.s3, amazon.s4, amazon.area)
amazon.mcwd.bar.2=figure5_barplot(amazon.annual.mcwd.sens.2[,,2], amazon.tmf.005, amazon.s1, amazon.s2, amazon.s3, amazon.s4, amazon.area)

amazon.tair.bar.1=figure5_barplot(amazon.annual.tair.sens.2[,,1], amazon.tmf.005, amazon.s1, amazon.s2, amazon.s3, amazon.s4, amazon.area)
amazon.tair.bar.2=figure5_barplot(amazon.annual.tair.sens.2[,,2], amazon.tmf.005, amazon.s1, amazon.s2, amazon.s3, amazon.s4, amazon.area)


africa_chirps_select_zscore_2016=africa.mcwd.hydro.zscore.2016 #africa_chirps_annual_mcwd_june_zscore_2016
#africa_chirps_select_zscore_2016[africa_chirps_select_zscore_2016>0]=NA

africa.s1=africa_chirps_select_zscore_2016[1:600, 1:400]>=(-1) & (africa.hydro.conform.2016<0)  # non severe dry 
africa.s2=africa_chirps_select_zscore_2016[1:600, 1:400]>=(-1) & (africa.hydro.conform.2016>0)  # non severe wet 
africa.s3=africa_chirps_select_zscore_2016[1:600, 1:400]<(-1) & (africa.hydro.conform.2016<0)  # severe dry 
africa.s4=africa_chirps_select_zscore_2016[1:600, 1:400]<(-1) & (africa.hydro.conform.2016>0)  # severe wet 

africa.mask5=array(NA, dim=c(600,400))
africa.mask5[africa.s1==1]=1
africa.mask5[africa.s2==1]=2
africa.mask5[africa.s3==1]=3
africa.mask5[africa.s4==1]=4
africa.mask5=africa.mask5*(africa.tmf.005.match>=0.95)
africa.mask5[africa.mask5==0]=NA

africa.mcwd.bar.1=figure5_barplot(congo.annual.mcwd.sens.2[,,1], africa.tmf.005.match, africa.s1, africa.s2, africa.s3, africa.s4, congo.area)
africa.mcwd.bar.2=figure5_barplot(congo.annual.mcwd.sens.2[,,2], africa.tmf.005.match, africa.s1, africa.s2, africa.s3, africa.s4, congo.area)

africa.tair.bar.1=figure5_barplot(congo.annual.tair.sens.2[,,1], africa.tmf.005.match, africa.s1, africa.s2, africa.s3, africa.s4, congo.area)
africa.tair.bar.2=figure5_barplot(congo.annual.tair.sens.2[,,2], africa.tmf.005.match, africa.s1, africa.s2, africa.s3, africa.s4, congo.area)



plot_bar_water=function(bar.1, bar.2, main.name, ylabel, ymin, ymax) {  # orange and red severe and nonsevere dry 
    cols=brewer.pal(4,'Spectral')[4:1]
    colos=c(cols[3], cols[2], cols[4], cols[1])
    colos[1]='#993404'
    colos[2]='#238443'
    colos[4]='#253494'
    par(mgp=c(3,1,0), mar=c(5,7,4,3))
    list1=c(bar.1, bar.2)
    plot(1:10,1:10, col='white', ylim=c(ymin, ymax), xlim=c(0.5,4.5), xlab='', ylab='', cex.lab=2, cex.axis=2,yaxs='i',las=1,bty='n', 
    main='', cex.main=2, xaxt='n')
    axis(1, at=1:4, labels=c(expression('R'[mild]^'dry'),expression('R'[mild]^'wet'),expression('R'[severe]^'dry'),expression('R'[severe]^'wet')),cex=2, cex.axis=2, lwd=0)
    mtext(side=2, text=ylabel, line=5, cex.axis=2, cex=1.5)
    mtext(side=3, text=main.name, cex=2, line=1)
    # xleft, ybottom, xright, ytop,
    rect(0.75, 0, 1, bar.1[1], border=colos[1], col='#c6dbef', lwd=5)
    rect(1, 0, 1.25, bar.2[1], border=colos[1], col='#4292c6', lwd=5)

    rect(1.75, 0, 2, bar.1[2], border=colos[2], col='#c6dbef', lwd=5)
    rect(2, 0, 2.25, bar.2[2], border=colos[2], col='#4292c6', lwd=5)

    rect(2.75, 0, 3, bar.1[3], border=colos[3], col='#c6dbef', lwd=5)
    rect(3, 0, 3.25, bar.2[3], border=colos[3], col='#4292c6', lwd=5)

    rect(3.75, 0, 4, bar.1[4], border=colos[4], col='#c6dbef', lwd=5)
    rect(4, 0, 4.25, bar.2[4], border=colos[4], col='#4292c6', lwd=5)
    abline(h=0)
}


plot_bar_tair=function(bar.1, bar.2, main.name, ylabel, ymin, ymax) {  # orange and red severe and nonsevere dry 
    cols=brewer.pal(4,'Spectral')[4:1]
    colos=c(cols[3], cols[2], cols[4], cols[1])
    colos[1]='#993404'
    colos[2]='#238443'
    colos[4]='#253494'
    par(mar=c(5,7,4,3))
    list1=c(bar.1, bar.2)
    plot(1:10,1:10, col='white', ylim=c(ymin, ymax), xlim=c(0.5,4.5), xlab='', ylab='', cex.lab=2, cex.axis=2,yaxs='i',las=1,bty='n',
    main='', cex.main=2, xaxt='n')
    #axis(1, at=1:4, labels=c('R1','R2','R3','R4'),cex=2, cex.axis=2)
    axis(1, at=1:4, labels=c(expression('R'[mild]^'dry'),expression('R'[mild]^'wet'),expression('R'[severe]^'dry'),expression('R'[severe]^'wet')),cex=2, cex.axis=2,lwd=0)
    mtext(side=2, text=ylabel, line=5, cex.axis=2, cex=1.5)
    mtext(side=3, text=main.name, cex=2, line=1)
    # xleft, ybottom, xright, ytop,
    rect(0.75, 0, 1, bar.1[1], border=colos[1], col='#fee391', lwd=5)
    rect(1, 0, 1.25, bar.2[1], border=colos[1], col='#fe9929', lwd=5)

    rect(1.75, 0, 2, bar.1[2], border=colos[2], col='#fee391', lwd=5)
    rect(2, 0, 2.25, bar.2[2], border=colos[2], col='#fe9929', lwd=5)

    rect(2.75, 0, 3, bar.1[3], border=colos[3], col='#fee391', lwd=5)
    rect(3, 0, 3.25, bar.2[3], border=colos[3], col='#fe9929', lwd=5)

    rect(3.75, 0, 4, bar.1[4], border=colos[4], col='#fee391', lwd=5)
    rect(4, 0, 4.25, bar.2[4], border=colos[4], col='#fe9929', lwd=5)
    abline(h=0)
}




ss1=(amazon.tmf.005>=0.95)*(amazon.mask5==1)*amazon.area
ss2=(amazon.tmf.005>=0.95)*(amazon.mask5==2)*amazon.area
ss3=(amazon.tmf.005>=0.95)*(amazon.mask5==3)*amazon.area
ss4=(amazon.tmf.005>=0.95)*(amazon.mask5==4)*amazon.area

total.area=(amazon.tmf.005>=0.95)*(!is.na(amazon.mask5))*amazon.area
total.area=(amazon.tmf.005>=0.95)*amazon.area

amazon.s1.frac=sum(ss1, na.rm=TRUE)/sum(total.area, na.rm=TRUE)*100
amazon.s2.frac=sum(ss2, na.rm=TRUE)/sum(total.area, na.rm=TRUE)*100
amazon.s3.frac=sum(ss3, na.rm=TRUE)/sum(total.area, na.rm=TRUE)*100
amazon.s4.frac=sum(ss4, na.rm=TRUE)/sum(total.area, na.rm=TRUE)*100


ss1=(africa.tmf.005.match>=0.95)*(africa.mask5==1)*congo.area
ss2=(africa.tmf.005.match>=0.95)*(africa.mask5==2)*congo.area
ss3=(africa.tmf.005.match>=0.95)*(africa.mask5==3)*congo.area
ss4=(africa.tmf.005.match>=0.95)*(africa.mask5==4)*congo.area

total.area=(africa.tmf.005.match>=0.95)*(!is.na(africa.mask5))*congo.area
total.area=(africa.tmf.005.match>=0.95)*congo.area

africa.s1.frac=sum(ss1, na.rm=TRUE)/sum(total.area, na.rm=TRUE)*100
africa.s2.frac=sum(ss2, na.rm=TRUE)/sum(total.area, na.rm=TRUE)*100
africa.s3.frac=sum(ss3, na.rm=TRUE)/sum(total.area, na.rm=TRUE)*100
africa.s4.frac=sum(ss4, na.rm=TRUE)/sum(total.area, na.rm=TRUE)*100

library(RColorBrewer)
library(maps)
library(fields)

tiff('figure_5barplot_tair_0112_hydro_oct.tiff', unit='in',width=14, height=7, res=200)
par(mfrow=c(2,3))

plot_bar_water(amazon.mcwd.bar.1*(-1), amazon.mcwd.bar.2*(-1), expression('(a)'~gamma[MCWD]), 'Sensitivity (dB/100mm)', -0.045, 0.022)
text(2.15, 0.015, '1992-2006', cex=2)
text(3.85, 0.015, '2006-2020', cex=2)
arrows(2.875,0,2.45,0.01, length=0.05)
arrows(3.125,0,3.55,0.01, length=0.05)

plot_bar_tair(amazon.tair.bar.1, amazon.tair.bar.2, expression('(b)'~gamma[Tair]), 'Sensitivity (dB/°C)', -0.12, 0.02)
#plot_bar(amazon.tair.bar.1, amazon.tair.bar.2, expression('(c)'~gamma[VPD]), 'Sensitivity (dB/hPa)')


#par(fig=c(0.15, 0.25, 0.6, 0.74),mar=c(0,0,0,0),new=TRUE)
par(mar=c(5,5,4,7))
#colos=brewer.pal(4,'Spectral')
cols=brewer.pal(4,'Spectral')[4:1]
colos=c(cols[3], cols[2], cols[4], cols[1])
colos[1]='#993404'
colos[2]='#238443'
colos[4]='#253494'
image(seq(-80,-50,by=1/20), seq(-20,10,by=1/20), amazon.mask5, xlab='', ylab='',  col=colos, main='', cex.main=2, cex.axis=2, cex.lab=2, las=1)
mtext(side=3, '(c) Amazon', cex=1.9, line=1)
map(database = "world",
        mar=c(1,1,1,1),col='gray', panel.first = grid(),add=TRUE,lwd=2)

#legend('bottom', legend=c(paste0('R1\n',round(amazon.s1.frac),'%'),
#paste0('R2\n', round(amazon.s2.frac),'%'),
#paste0('R3\n', round(amazon.s3.frac),'%'),
#paste0('R4\n', round(amazon.s4.frac),'%')), 
#text.col=colos, cex=3, bty='n', horiz=TRUE)

legend('bottomleft', legend=c(paste0(round(amazon.s1.frac),'%'),
paste0( round(amazon.s2.frac),'%'),
paste0( round(amazon.s3.frac),'%'),
paste0( round(amazon.s4.frac),'%')), 
text.col=colos, cex=2, bty='n', ncol=2)

plot_bar_water(africa.mcwd.bar.1*(-1), africa.mcwd.bar.2*(-1), expression('(d)'~gamma[MCWD]) , 'Sensitivity (dB/100mm)', -0.045, 0.022)
cols=brewer.pal(4,'Spectral')[4:1]
colos=c(cols[3], cols[2], cols[4], cols[1])
colos[1]='#993404'
#colos[1]='#f781bf'
colos[2]='#238443'
colos[4]='#253494'

plot_bar_tair(africa.tair.bar.1, africa.tair.bar.2, expression('(e)'~gamma[Tair]), 'Sensitivity (dB/°C)', -0.12, 0.02)
#plot_bar(africa.vpd.bar.1, africa.vpd.bar.2, expression('(f)'~gamma[VPD]), 'Sensitivity (dB/hPa)')
legend('bottomleft', y.intersp=1.2,
legend=c(expression('R'[mild]^'dry'~': Dry season (Z>-1)'),expression('R'[mild]^'wet'~': Wet season (Z>-1)'),
expression('R'[severe]^'dry'~': Dry season (Z<-1)'), expression('R'[severe]^'wet'~': Wet season (Z<-1)')), text.col=colos, cex=2, bty='n')
 

par(mar=c(5,5,4,7))
#par(fig=c(0.35, 0.45, 0.12, 0.22),mar=c(0,0,0,0),new=TRUE)
image(seq(5,35,by=1/20), seq(-10,10,by=1/20), africa.mask5, xlab='', ylab='',   col=colos, main='', cex.main=2, cex.axis=2, cex.lab=2, las=1)
mtext(side=3, '(f) C. Africa', cex=1.9, line=1)
map(database = "world",
        mar=c(1,1,1,1),col='gray', panel.first = grid(),add=TRUE,lwd=2)

#legend('bottom', legend=c(paste0('R1\n',round(africa.s1.frac),'%'),
#paste0('R2\n', round(africa.s2.frac),'%'),
#paste0('R3\n', round(africa.s3.frac),'%'),
#paste0('R4\n', round(africa.s4.frac),'%')), 
#text.col=colos, cex=3, bty='n', horiz=TRUE)

legend('bottomleft', legend=c(paste0(round(africa.s1.frac),'%'),
paste0( round(africa.s2.frac),'%'),
paste0( round(africa.s3.frac),'%'),
paste0( round(africa.s4.frac),'%')), 
text.col=colos, cex=2, bty='n',  ncol=2)

dev.off()

