
load('SST_nino34.Rdata')
load('/users/yyao/Documents/futureAmazon/area_amazon_congo_005.Rdata')
load(file='/users/yyao/Documents/futureAmazon/amazon_fraction_tmf_deg01.Rdata')
load(file='/users/yyao/Documents/futureAmazon/amazon_fraction_tmf_deg005.Rdata')

#load(file='/users/yyao/Documents/futureAmazon/africa_fraction_tmf_deg05.Rdata')
load(file='/users/yyao/Documents/futureAmazon/africa_fraction_tmf_deg01_June22.Rdata')
load(file='/users/yyao/Documents/futureAmazon/africa_fraction_tmf_deg005_June22.Rdata')

# 10 to 35 -10 to 10 
#africa.tmf.05.match=africa.tmf.05[1:50, ]
africa.tmf.01.match=africa.tmf.01[51:350, ]
africa.tmf.005.match=africa.tmf.005[101:700, ]  # start from 5 

load(file='/users/yyao/Documents/futureAmazon/africa_chirps_hydro_month_mcwd_zscore_newbase1992-2020_0722.Rdata')
africa.1983.month.mcwd=array(NA, dim=c(600,400,24))
africa.1998.month.mcwd=array(NA, dim=c(600,400,24))
africa.2016.month.mcwd=array(NA, dim=c(600,400,24))
africa.2024.month.mcwd=array(NA, dim=c(600,400,24))

africa.1983.month.mcwd[, , 1:9]=africa.chirps.month.mcwd.zscore[1:600, 1:400, 1982-1982+1, 4:12]
africa.1983.month.mcwd[, , 10:21]=africa.chirps.month.mcwd.zscore[1:600, 1:400, 1983-1982+1, 1:12]
africa.1983.month.mcwd[, , 22:24]=africa.chirps.month.mcwd.zscore[1:600, 1:400, 1984-1982+1, 1:3]


africa.1998.month.mcwd[, , 1:9]=africa.chirps.month.mcwd.zscore[1:600, 1:400, 1997-1982+1, 4:12]
africa.1998.month.mcwd[, , 10:21]=africa.chirps.month.mcwd.zscore[1:600, 1:400, 1998-1982+1, 1:12]
africa.1998.month.mcwd[, , 22:24]=africa.chirps.month.mcwd.zscore[1:600, 1:400, 1999-1982+1, 1:3]


africa.2016.month.mcwd[, , 1:9]=africa.chirps.month.mcwd.zscore[1:600, 1:400, 2015-1982+1, 4:12]
africa.2016.month.mcwd[, , 10:21]=africa.chirps.month.mcwd.zscore[1:600, 1:400, 2016-1982+1, 1:12]
africa.2016.month.mcwd[, , 22:24]=africa.chirps.month.mcwd.zscore[1:600, 1:400, 2017-1982+1, 1:3]


africa.2024.month.mcwd[, , 1:9]=africa.chirps.month.mcwd.zscore[1:600, 1:400, 2023-1982+1, 4:12]
africa.2024.month.mcwd[, , 10:21]=africa.chirps.month.mcwd.zscore[1:600, 1:400, 2024-1982+1, 1:12]
#africa.2024.month.mcwd[, , 22:24]=africa.chirps.month.mcwd.zscore[1:600, 101:500, 2024-1982+1, 1:3]

#### running average 

africa.area.1983=array(NA, dim=c(4, 24))
africa.area.1998=array(NA, dim=c(4, 24))
africa.area.2016=array(NA, dim=c(4, 24))
africa.area.2024=array(NA, dim=c(4, 24))
year=1982
for (month in 1:24) {
    zscore=africa.1983.month.mcwd[, , month]*(africa.tmf.005.match>=0.95)
    zscore.mask=((zscore>(-1.645)) & (zscore<=(-1)))*congo.area
    africa.area.1983[1, month]=sum(zscore.mask, na.rm=TRUE)
    zscore.mask=((zscore>(-1.96)) & (zscore<=(-1.645)))*congo.area
    africa.area.1983[2, month]=sum(zscore.mask, na.rm=TRUE)
    zscore.mask=((zscore>(-2.576)) & (zscore<=(-1.96)))*congo.area
    africa.area.1983[3, month]=sum(zscore.mask, na.rm=TRUE)
    zscore.mask=(zscore<=(-2.576))*congo.area
    africa.area.1983[4, month]=sum(zscore.mask, na.rm=TRUE)
}

year=1997
for (month in 1:24) {
    zscore=africa.1998.month.mcwd[, , month]*(africa.tmf.005.match>=0.95)
    zscore.mask=((zscore>(-1.645)) & (zscore<=(-1)))*congo.area
    africa.area.1998[1, month]=sum(zscore.mask, na.rm=TRUE)
    zscore.mask=((zscore>(-1.96)) & (zscore<=(-1.645)))*congo.area
    africa.area.1998[2, month]=sum(zscore.mask, na.rm=TRUE)
    zscore.mask=((zscore>(-2.576)) & (zscore<=(-1.96)))*congo.area
    africa.area.1998[3, month]=sum(zscore.mask, na.rm=TRUE)
    zscore.mask=(zscore<=(-2.576))*congo.area
    africa.area.1998[4, month]=sum(zscore.mask, na.rm=TRUE)
}

year=2015
for (month in 1:24) {
    zscore=africa.2016.month.mcwd[, , month]*(africa.tmf.005.match>=0.95)
    zscore.mask=((zscore>(-1.645)) & (zscore<=(-1)))*congo.area
    africa.area.2016[1, month]=sum(zscore.mask, na.rm=TRUE)
    zscore.mask=((zscore>(-1.96)) & (zscore<=(-1.645)))*congo.area
    africa.area.2016[2, month]=sum(zscore.mask, na.rm=TRUE)
    zscore.mask=((zscore>(-2.576)) & (zscore<=(-1.96)))*congo.area
    africa.area.2016[3, month]=sum(zscore.mask, na.rm=TRUE)
    zscore.mask=(zscore<=(-2.576))*congo.area
    africa.area.2016[4, month]=sum(zscore.mask, na.rm=TRUE)
}

for (month in 1:24) {
    zscore=africa.2024.month.mcwd[, , month]*(africa.tmf.005.match>=0.95)
    zscore.mask=((zscore>(-1.645)) & (zscore<=(-1)))*congo.area
    africa.area.2024[1, month]=sum(zscore.mask, na.rm=TRUE)
    zscore.mask=((zscore>(-1.96)) & (zscore<=(-1.645)))*congo.area
    africa.area.2024[2, month]=sum(zscore.mask, na.rm=TRUE)
    zscore.mask=((zscore>(-2.576)) & (zscore<=(-1.96)))*congo.area
    africa.area.2024[3, month]=sum(zscore.mask, na.rm=TRUE)
    zscore.mask=(zscore<=(-2.576))*congo.area
    africa.area.2024[4, month]=sum(zscore.mask, na.rm=TRUE)
}


load(file='/users/yyao/Documents/futureAmazon/amazon_chirps_hydro_month_mcwd_zscore_600_600_newbase1992-2020_0722.Rdata')

amazon.1983.month.mcwd=array(NA, dim=c(600,600,24))
amazon.1998.month.mcwd=array(NA, dim=c(600,600,24))
amazon.2016.month.mcwd=array(NA, dim=c(600,600,24))
amazon.2024.month.mcwd=array(NA, dim=c(600,600,24))

amazon.1983.month.mcwd[, , 1:9]=amazon.chirps.month.mcwd.zscore[, , 1982-1982+1, 4:12]
amazon.1983.month.mcwd[, , 10:21]=amazon.chirps.month.mcwd.zscore[, , 1983-1982+1, 1:12]
amazon.1983.month.mcwd[, , 22:24]=amazon.chirps.month.mcwd.zscore[, , 1984-1982+1, 1:3]


amazon.1998.month.mcwd[, , 1:9]=amazon.chirps.month.mcwd.zscore[, , 1997-1982+1, 4:12]
amazon.1998.month.mcwd[, , 10:21]=amazon.chirps.month.mcwd.zscore[, , 1998-1982+1, 1:12]
amazon.1998.month.mcwd[, , 22:24]=amazon.chirps.month.mcwd.zscore[, , 1999-1982+1, 1:3]


amazon.2016.month.mcwd[, , 1:9]=amazon.chirps.month.mcwd.zscore[, , 2015-1982+1, 4:12]
amazon.2016.month.mcwd[, , 10:21]=amazon.chirps.month.mcwd.zscore[, , 2016-1982+1, 1:12]
amazon.2016.month.mcwd[, , 22:24]=amazon.chirps.month.mcwd.zscore[, , 2017-1982+1, 1:3]


amazon.2024.month.mcwd[, , 1:9]=amazon.chirps.month.mcwd.zscore[, , 2023-1982+1, 4:12]
amazon.2024.month.mcwd[, , 10:21]=amazon.chirps.month.mcwd.zscore[, , 2024-1982+1, 1:12]
#amazon.2024.month.mcwd[, , 22:24]=amazon.chirps.month.mcwd.zscore[, , 2017-1982+1, 1:3]


amazon.area.1983=array(NA, dim=c(4, 24))
amazon.area.1998=array(NA, dim=c(4, 24))
amazon.area.2016=array(NA, dim=c(4, 24))
amazon.area.2024=array(NA, dim=c(4, 24))

year=1982
for (month in 1:24) {
    zscore=amazon.1983.month.mcwd[, , month]*(amazon.tmf.005>=0.95)
    zscore.mask=((zscore>(-1.645)) & (zscore<=(-1)))*amazon.area
    amazon.area.1983[1, month]=sum(zscore.mask, na.rm=TRUE)
    zscore.mask=((zscore>(-1.96)) & (zscore<=(-1.645)))*amazon.area
    amazon.area.1983[2, month]=sum(zscore.mask, na.rm=TRUE)
    zscore.mask=((zscore>(-2.576)) & (zscore<=(-1.96)))*amazon.area
    amazon.area.1983[3, month]=sum(zscore.mask, na.rm=TRUE)
    zscore.mask=(zscore<=(-2.576))*amazon.area
    amazon.area.1983[4, month]=sum(zscore.mask, na.rm=TRUE)
}

year=1997
for (month in 1:24) {
    zscore=amazon.1998.month.mcwd[, , month]*(amazon.tmf.005>=0.95)
    zscore.mask=((zscore>(-1.645)) & (zscore<=(-1)))*amazon.area
    amazon.area.1998[1, month]=sum(zscore.mask, na.rm=TRUE)
    zscore.mask=((zscore>(-1.96)) & (zscore<=(-1.645)))*amazon.area
    amazon.area.1998[2, month]=sum(zscore.mask, na.rm=TRUE)
    zscore.mask=((zscore>(-2.576)) & (zscore<=(-1.96)))*amazon.area
    amazon.area.1998[3, month]=sum(zscore.mask, na.rm=TRUE)
    zscore.mask=(zscore<=(-2.576))*amazon.area
    amazon.area.1998[4, month]=sum(zscore.mask, na.rm=TRUE)
}

year=2015
for (month in 1:24) {
    zscore=amazon.2016.month.mcwd[, , month]*(amazon.tmf.005>=0.95)
    zscore.mask=((zscore>(-1.645)) & (zscore<=(-1)))*amazon.area
    amazon.area.2016[1, month]=sum(zscore.mask, na.rm=TRUE)
    zscore.mask=((zscore>(-1.96)) & (zscore<=(-1.645)))*amazon.area
    amazon.area.2016[2, month]=sum(zscore.mask, na.rm=TRUE)
    zscore.mask=((zscore>(-2.576)) & (zscore<=(-1.96)))*amazon.area
    amazon.area.2016[3, month]=sum(zscore.mask, na.rm=TRUE)
    zscore.mask=(zscore<=(-2.576))*amazon.area
    amazon.area.2016[4, month]=sum(zscore.mask, na.rm=TRUE)
}

year=2024
for (month in 1:24) {
    zscore=amazon.2024.month.mcwd[, , month]*(amazon.tmf.005>=0.95)
    zscore.mask=((zscore>(-1.645)) & (zscore<=(-1)))*amazon.area
    amazon.area.2024[1, month]=sum(zscore.mask, na.rm=TRUE)
    zscore.mask=((zscore>(-1.96)) & (zscore<=(-1.645)))*amazon.area
    amazon.area.2024[2, month]=sum(zscore.mask, na.rm=TRUE)
    zscore.mask=((zscore>(-2.576)) & (zscore<=(-1.96)))*amazon.area
    amazon.area.2024[3, month]=sum(zscore.mask, na.rm=TRUE)
    zscore.mask=(zscore<=(-2.576))*amazon.area
    amazon.area.2024[4, month]=sum(zscore.mask, na.rm=TRUE)
}


amazon.area.1983[is.na(amazon.area.1983)]=0
amazon.area.1998[is.na(amazon.area.1998)]=0
amazon.area.2016[is.na(amazon.area.2016)]=0
amazon.area.2024[is.na(amazon.area.2024)]=0

africa.area.1983[is.na(africa.area.1983)]=0
africa.area.1998[is.na(africa.area.1998)]=0
africa.area.2016[is.na(africa.area.2016)]=0
africa.area.2024[is.na(africa.area.2024)]=0

amazon.area.1983.frac=amazon.area.1983/sum(amazon.area*(amazon.tmf.005>=0.95), na.rm=TRUE)*100
amazon.area.1998.frac=amazon.area.1998/sum(amazon.area*(amazon.tmf.005>=0.95), na.rm=TRUE)*100
amazon.area.2016.frac=amazon.area.2016/sum(amazon.area*(amazon.tmf.005>=0.95), na.rm=TRUE)*100
amazon.area.2024.frac=amazon.area.2024/sum(amazon.area*(amazon.tmf.005>=0.95), na.rm=TRUE)*100

africa.area.1983.frac=africa.area.1983/sum(congo.area*(africa.tmf.005.match>=0.95), na.rm=TRUE)*100
africa.area.1998.frac=africa.area.1998/sum(congo.area*(africa.tmf.005.match>=0.95), na.rm=TRUE)*100
africa.area.2016.frac=africa.area.2016/sum(congo.area*(africa.tmf.005.match>=0.95), na.rm=TRUE)*100
africa.area.2024.frac=africa.area.2024/sum(congo.area*(africa.tmf.005.match>=0.95), na.rm=TRUE)*100




plot_area=function(area.1983, ymax, title) {
    plot(area.1983[4, ], ylim=c(0,ymax), xaxt='n', yaxt='n',xlab='', ylab='Area fraction (%)', cex.axis=2, cex.lab=2, yaxs='i', main=title, cex.main=2, las=1)
    xx=1:24
    Con.low=rep(0,24)
    Con.high=area.1983[4, ]
    polygon(c(xx,rev(xx)),c(Con.low,rev(Con.high)),col=colos[4],border=NA)
    Con.low=area.1983[4, ]
    Con.high=area.1983[4, ]+area.1983[3, ]
    polygon(c(xx,rev(xx)),c(Con.low,rev(Con.high)),col=colos[3],border=NA)
    Con.low=area.1983[4, ]+area.1983[3, ]
    Con.high=area.1983[4, ]+area.1983[3, ]+area.1983[2, ]
    polygon(c(xx,rev(xx)),c(Con.low,rev(Con.high)),col=colos[2],border=NA)
    Con.low=area.1983[4, ]+area.1983[3, ]+area.1983[2, ]
    Con.high=area.1983[4, ]+area.1983[3, ]+area.1983[2, ]+area.1983[1, ]
    polygon(c(xx,rev(xx)),c(Con.low,rev(Con.high)),col=colos[1],border=NA)
    #axis(1, at=1:24, labels=c('J','','','A','','','J','','','O','','','J','','','A','','','J','','','O','',''), cex.axis=2, cex=2)
    axis(1, at=seq(1,24,by=3), labels=c('J','A','J','O','J','A','J','O'), cex.axis=2, cex=2)
    abline(v=12, lwd=1, lty=1)
    axis(2, at=seq(0,ymax,by=20), labels=seq(0,ymax,by=20), las=1, cex.axis=2, cex=2)

}


#oni=read.csv('oni.csv')
nino1996=nino$ANOM.3[nino$YR==1996]
nino1997=nino$ANOM.3[nino$YR==1997]
nino1998=nino$ANOM.3[nino$YR==1998]

nino2014=nino$ANOM.3[nino$YR==2014]
nino2015=nino$ANOM.3[nino$YR==2015]
nino2016=nino$ANOM.3[nino$YR==2016]


tiff('figure1_0112.tiff', width=1000, height=800)

colos=c('#fed976','#fd8d3c','#e31a1c','#800026')
par(mfrow=c(2, 2))

par(fig=c(0.05, 0.925, 0.7, 0.975), mar=c(2,5,2,2), new=TRUE)  
nino.index=c(nino2015, nino2016)
plot(1:24, nino.index, xlim=c(0.4,24.5), xlab='', ylab='SST anomalies (°C)', ylim=c(-2,3), col='white', xaxt='n', yaxt='n',las=1,
cex.axis=2, cex.lab=2, main='(a) Niño3.4 SST anomalies index', cex.main=2)
lines(1:24, nino.index, lwd=2, col='red')
points(1:24, nino.index, pch=1, cex=2, lwd=2, col='red')
#axis(side = 4, at = pretty(range(oni.index)), cex.axis=2)    
nino.index=c(nino1997, nino1998)
lines(1:24, nino.index, lwd=2, lty=1 )
points(1:24, nino.index, pch=3, cex=2, lwd=2)
abline(h=0, lty=2)
abline(h=1, lty=2)
abline(h=2, lty=2)
abline(h=-1, lty=2)

#lines(1:24, c(nino.mean, nino.mean), lwd=1, lty=1)
#points(1:24, c(nino.mean, nino.mean), pch=5, cex=2, lwd=2)

#abline(h=-1, lwd=1, lty=2, col='red')
#abline(h=1, lwd=1, lty=2, col='red')
abline(v=12, lwd=1, lty=1)
abline(v=6, lty=2, lwd=1)
abline(v=17, lty=2, lwd=1)
axis(1, at=1:24, labels=c('J','F','M','A','M','J','J','A','S','O','N','D','J','F','M','A','M','J','J','A','S','O','N','D'), cex.axis=2, cex=2)
axis(2, at=-2:3, labels=seq(-2,3,by=1),las=1, cex.axis=2, cex=2)
text(6.5, -1, 'June', cex=2, srt=-90)
text(17.5, -1, 'May', cex=2, srt=-90)
text(13,-1,'Dec', cex=2, srt=-90)
abline(h=29)
legend('topright', legend=c('1997/98','2015/16'), pch=c(3,1),col=c('black','red'), text.col=c('black','red'), bty='n', cex=2)

par(mar=c(5,5,3,5))
par(fig=c(0.05, 0.47, 0.4, 0.7), mar=c(2,5,2,0), new=TRUE) 
plot_area(amazon.area.1998.frac, 70, '(b) Amazon 1997/98')
#text(20,40, '1997/98', cex=2)

abline(v=6, lty=2, lwd=1)
abline(v=17, lty=2, lwd=1)
legend('topright',legend=c('Mild','Moderate','Severe','Extreme'),text.col=colos, bty='n', cex=2)


par(mar=c(5,5,3,5))
par(fig=c(0.5, 0.9, 0.4, 0.7), mar=c(2,2,2,0), new=TRUE) 
plot_area(amazon.area.2016.frac, 70, '(c) Amazon 2015/16')
#text(20,40, '2015/16', cex=2)
#text(13.5,41.25,'Dec', cex=2, srt=-90)
abline(v=6, lty=2, lwd=1)
abline(v=17, lty=2, lwd=1)

#text(7, 41.25, 'June', cex=2, srt=-90)
#text(18, 41.25, 'May', cex=2, srt=-90)


par(mar=c(5,5,3,5))
par(fig=c(0.05, 0.47, 0.05, 0.35), mar=c(2,5,2,0), new=TRUE) 
plot_area(africa.area.1998.frac, 70, '(d) C. Africa 1997/98')
#text(20,64, '1997/98', cex=2)
#text(13.5, 66, 'Dec', cex=2, srt=-90)
abline(v=6, lty=2, lwd=1)
abline(v=17, lty=2, lwd=1)
#text(7, 66, 'June', cex=2, srt=-90)
#text(18, 66, 'May', cex=2, srt=-90)


par(mar=c(5,5,3,5))
par(fig=c(0.5, 0.9, 0.05, 0.35), mar=c(2,2,2,0), new=TRUE) 
plot_area(africa.area.2016.frac, 70, '(e) C. Africa 2015/16')
#text(20,64, '2015/16', cex=2)
#text(13.5,66,'Dec', cex=2, srt=-90)
abline(v=6, lty=2, lwd=1)
abline(v=17, lty=2, lwd=1)
#text(7, 66, 'June', cex=2, srt=-90)
#text(18, 66, 'May', cex=2, srt=-90)

dev.off()


