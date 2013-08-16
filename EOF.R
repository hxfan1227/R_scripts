################### EOF ANALYSIS ####################
##Concatenates all lon sections of annual/winter/summer means and calculates EOFs 

### parameters
set<-	3	##if dataset is too large, this breaks it into n sets along lon. 
##have as many 'data' objects as there are sets!
data1 <- "CCSM4_pr_piC_clim1.nc"	##netCDF filename to extract data from
data2 <- "CCSM4_pr_piC_clim2.nc"	##netCDF filename to extract data from
data3 <- "CCSM4_pr_piC_clim3.nc"	##netCDF filename to extract data from

orig<-"pr_Amon_CCSM4_piControl_r1i1p1_025001-050012.nc" ##original dataset; used to define dimensions of new netCDF

var<- "annual_mean"		##variable to extract and create, suggested format "xxx_mean" where xxx are months/seasons 
units<-"mm/day" 		##variable units
longname<-"Annual Mean Precipitation 250-1300"

start<-250			##dataset starting year
end<- 1300			##dataset ending year
yy<-length(start:end)		##dataset time length 

eofnc<-	"CCSM4_pr_piC_annualmeanEOF.nc"		##new netCDF for EOF analysis
eof1long<-"EOF1 250-1300"		##new variable long name, suggested format 
eof2long<-"EOF2 250-1300"		##new variable long name, suggested format
eof3long<-"EOF3 250-1300"		##new variable long name, suggested format
eoffile<- "CCSM4_pr_piC_EOF.csv"	##where to store EOF/PC output 
	## use"[model]_[variable]_[run]_[meantype]EOF.csv"
varfile<-"CCSM4_pr_piC_var.csv"		##where to store var/corr output


######## Create new netCDF file for full means and EOFs 
##define dimensions; lon from original dataset to get full range
nc<-open.ncdf(orig, write=T)
lon<-nc$dim$lon
close.ncdf(nc)
##get lat and time from calculated datasets (should be same for all lon sets)
nc<-open.ncdf(data1, write=T)
lat<-nc$dim$lat
time<-nc$dim$time
time$vals<-1:yy


close.ncdf(nc)

##define variables
meanvar<-var.def.ncdf(var,units,dim=list(lon,lat,time),missval=NA,longname=longname)

##create netCDF file
nc<-create.ncdf(eofnc,meanvar)


close.ncdf(nc)

###2. ####### Concatenate means to get full lon 
##get lon1 set
nc<-open.ncdf(data1, write=T)
nc1<-get.var.ncdf(nc,var,start=c(1,1,1),count=c(-1,-1,yy)) 
lon1<-nc$dim$lon1$len
close.ncdf(nc)
##get lon2 set
nc<-open.ncdf(data2, write=T)
nc2<-get.var.ncdf(nc,var,start=c(1,1,1),count=c(-1,-1,yy)) 
lon2<-nc$dim$lon2$len
close.ncdf(nc)
##get lon3 set
nc<-open.ncdf(data3, write=T)
nc3<-get.var.ncdf(nc,var,start=c(1,1,1),count=c(-1,-1,yy)) 
lon3<-nc$dim$lon3$len
close.ncdf(nc)
##combine all lon sets
ncmean<-abind(nc1,nc2,nc3,along=1)  

##write means to netCDF
nc<-open.ncdf(eofnc,write=T)
put.var.ncdf(nc,var,ncmean)

close.ncdf(nc)


###3. ########### Calculate EOFs
##get mean dataset
nc<-open.ncdf(eofnc,write=T)
ncmean<-get.var.ncdf(nc,var)
lon<-nc$dim$lon$len
lat<-nc$dim$lat$len

close.ncdf(nc)
##detrend data
nc<-wrap(ncmean,map=list(NA,3)) #make a matrix [lon*lat,time]
nc<-detrend(nc) #detrend 
nc<-t(nc) #matrix [t,lon*lat] for calculating EOF

##for p>>n EOF and PC
L<-nc%*%t(nc) #small covariance matrix [t,t]
eig<-eigen(L) #eig$vectors and eig$values
lam<-diag(eig$values, nrow = length(eig$values)) #diagonal matrix Lambda (=values); equal to those of R=t(nc)%*%nc
B<-eig$vectors #eigenvectors B
D<-t(nc)%*%B #proportional EOFs
eof<-sweep(D, 2, sqrt(colSums(D^2)), FUN="/") #normalise EOFs

#write first 10 EOFs to .csv eoffile
write.table(eof[,1:10],eoffile,col.names=list("EOF1", "EOF2","EOF3","EOF4","EOF5","EOF6","EOF7","EOF8","EOF9","EOF10"))

##EOF1
eof1<-eof[,1] #1st EOF
pc1<-nc%*%eof1 #PC1
corr1<-cor(pc1,nc,method="spearman") #correlation1
dim(corr1)<-c(lon,lat) 
varexp1<-(corr1^2)*100

##EOF2
eof2<-eof[,2] #2nd EOF
pc2<-nc%*%eof2 #PC2
corr2<-cor(pc2,nc,method="spearman") #correlation2
dim(corr2)<-c(lon,lat) 
varexp2<-(corr2^2)*100

##EOF3
eof3<-eof[,3] #3rd EOF
pc3<-nc%*%eof3 #PC3
corr3<-cor(pc3,nc,method="spearman") #correlation3
dim(corr3)<-c(lon,lat) 
varexp3<-(corr3^2)*100

##combine PCs into matrix
PC<-abind(pc1,pc2,pc3,along=2)
#write PC to .csv eoffile
write.table(PC,eoffile,append=T,col.names=list("PC1","PC2","PC3"))

##variance explained by each EOF
variance<-diag(lam)/tr(lam) 
variance[1:10]*100 #variance % for first 10 EOFs

#write variance to .csv varfile
write.table(variance[1:10],varfile,col.names="Variance")

##combine correlations into matrix
correlation<-abind(corr1,corr2,corr3,along=2)
#define colnames (latitudes)
cols<-gl(3,lat,labels=list("Correlation1","Correlation2","Correlation3"))

#write correlation to .csv varfile
write.table(correlation,varfile,append=T,col.names=cols)





























################# is this necessary??? ################
#switch lon
east<-round(lon/2)
west<-lon-east
nc.e<-nc.1[,,1:east]
nc.w<-nc.1[,,west:lon]
nc1<-abind(nc.w,nc.e,along=3)

rm(nc.1,nc.e,nc.w,east,west)
#######################################################


############MAPS
##change EOF dimensions
dim(eof.1k)<-c(192,288) #1st EOF, in [lat,lon] last millennium
eof.1k<-t(eof.1k) #change to [lon,lat]

dim(eof.2k)<-c(192,288) #1st EOF, in [lat,lon] last millennium
eof.2k<-t(eof.2k) #change to [lon,lat]


#map settings
ccsm<-open.ncdf("pr_Amon_CCSM4_midHolocene_r1i1p1_100001-130012.nc") #extract lat, lon vectors
lon<-get.var.ncdf(ccsm,"lon") 
lat<-get.var.ncdf(ccsm,"lat")

data(world2HiresMapEnv) #map contours
colour<-colorRampPalette(c("red","white","blue"), space="rgb")
colour2<-colorRampPalette(c("white","yellow","orange","red"), space="rgb")



##past millennium
#global EOF1
png(filename="EOF_1_1k.png", width=930,height=560) ##EOF1
filled.contour(lon, lat, eof.1k*100, nlevels=40, zlim=c(-4,4), color.palette=colour, plot.title=title(main="1st EOF Last Millennium"),plot.axes={maps::map(database="world2Hires",interior=T,add=T)})
dev.off()

png(filename="Correlation_1_1k.png", width=930,height=560) ## correlation1
filled.contour(lon, lat, corr.1k, nlevels=40, zlim=c(-1,1), color.palette=colour, plot.title=title(main="1st EOF Correlation Last Millennium"),plot.axes={maps::map(database="world2Hires",interior=T,add=T)})
dev.off()

png(filename="Variance_1_1k.png", width=930,height=560) ## variance explained1
filled.contour(lon, lat, varexp.1k, nlevels=40, zlim=c(0,100), color.palette=colour2, plot.title=title(main="Variance Explained by EOF1 Last Millennium"),plot.axes={maps::map(database="world2Hires",interior=T,add=T)})
dev.off()

png("PC1_Var_1k.png", width=930,height=560) # PC1 and variance
plot.new()
par(mfrow=c(2,1)) #2 figures in 2 rows, 1 column
plot.ts(pc.1k,xlab="Year AD",ylab="",main="PC1 Last Millennium")
plot(var.1k[1:10]*100, main="Variances",xlab="EOF",ylab="Variance (%)",type="l")
dev.off()

#global EOF2
png(filename="EOF_2_1k.png", width=930,height=560) ##EOF1
filled.contour(lon, lat, eof.2k*100, nlevels=40, zlim=c(-4,4), color.palette=colour, plot.title=title(main="2nd EOF Last Millennium"),plot.axes={maps::map(database="world2Hires",interior=T,add=T)})
dev.off()

png(filename="Correlation_2_1k.png", width=930,height=560) ## correlation1
filled.contour(lon, lat, corr.2k, nlevels=40, zlim=c(-1,1), color.palette=colour, plot.title=title(main="2nd EOF Correlation Last Millennium"),plot.axes={maps::map(database="world2Hires",interior=T,add=T)})
dev.off()

png(filename="Variance_2_1k.png", width=930,height=560) ## variance explained1
filled.contour(lon, lat, varexp.2k, nlevels=40, zlim=c(0,100), color.palette=colour2, plot.title=title(main="Variance Explained by EOF2 Last Millennium"),plot.axes={maps::map(database="world2Hires",interior=T,add=T)})948
dev.off()

png("PC2_Var_1k.png", width=930,height=560) # PC1 and variance
plot.new()
par(mfrow=c(2,1)) #2 figures in 2 rows, 1 column
plot.ts(pc.2k,xlab="Year AD",ylab="",main="PC2 Last Millennium")
plot(var.1k[1:10]*100, main="Variances",xlab="EOF",ylab="Variance (%)",type="l")
dev.off()


png("PC1_MCA_LIA_1k.png", width=930,height=560) # PC1 and variance
plot.new()
par(mfrow=c(2,1)) #2 figures in 2 rows, 1 column
plot.ts(pc.1k[150:450],xlim=c(0,450),ylim=c(-150,150),xlab="Year since AD 1000",ylab="",main="PC1 MCA")
plot.ts(pc.1k[550:1000],xlim=c(0,450),ylim=c(-150,150),xlab="Year since AD 1400",ylab="",main="PC1 LIA")
dev.off()

























