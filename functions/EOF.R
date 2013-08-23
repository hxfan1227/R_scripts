################### EOF ANALYSIS ####################
##Concatenates all lon sections of annual/winter/summer means and calculates EOFs 

### parameters
set<-	3	##if dataset is too large, this breaks it into n sets along lon. 
##have as many 'data' objects as there are sets!
data1 <- "CCSM4_pr_mh_clim1.nc"	##netCDF filename to extract data from
data2 <- "CCSM4_pr_mh_clim2.nc"	##netCDF filename to extract data from
data3 <- "CCSM4_pr_mh_clim3.nc"	##netCDF filename to extract data from

orig<-"pr_Amon_CCSM4_midHolocene_r1i1p1_100001-130012.nc" ##original dataset; used to define dimensions of new netCDF

varmean<- "annual_mean"		##variable to extract and create, suggested format "xxx_mean" where xxx are months/seasons 
vareof<-"eof"

units<-"mm/day" 		##variable units
longmean<-"Annual Mean Precipitation Mid Holocene"
longeof<-"EOF"

start<-1000			##dataset starting year
end<- 1300			##dataset ending year
yy<-length(start:end)		##dataset time length 

eofnc<-	"CCSM4_pr_mh_meanEOF.nc"		##new netCDF for EOF analysis
eof1long<-"EOF1 Mid Holocene"		##new variable long name, suggested format 
eof2long<-"EOF2 Mid Holocene"		##new variable long name, suggested format
eof3long<-"EOF3 Mid Holocene"		##new variable long name, suggested format
eoffile<- "CCSM4_pr_mh_EOF.csv"	##where to store EOF/PC output 
	## use"[model]_[variable]_[run]_[meantype]EOF.csv"
varfile<-"CCSM4_pr_mh_var.csv"		##where to store var/corr output


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
meanvar<-var.def.ncdf(varmean,units,dim=list(lon,lat,time),missval=NA,longname=longmean)
eofvar<-var.def.ncdf(vareof,units,dim=list(lon,lat,time),missval=NA,longname=longeof)

##create netCDF file
nc<-create.ncdf(eofnc,list(meanvar,eofvar))


close.ncdf(nc)

###2. ####### Concatenate means to get full lon 
##get lon1 set
nc<-open.ncdf(data1, write=T)
nc1<-get.var.ncdf(nc,varmean,start=c(1,1,1),count=c(-1,-1,yy)) 
lon1<-nc$dim$lon1$len
close.ncdf(nc)
##get lon2 set
nc<-open.ncdf(data2, write=T)
nc2<-get.var.ncdf(nc,varmean,start=c(1,1,1),count=c(-1,-1,yy)) 
lon2<-nc$dim$lon2$len
close.ncdf(nc)
##get lon3 set
nc<-open.ncdf(data3, write=T)
nc3<-get.var.ncdf(nc,varmean,start=c(1,1,1),count=c(-1,-1,yy)) 
lon3<-nc$dim$lon3$len
close.ncdf(nc)
##combine all lon sets
ncmean<-abind(nc1,nc2,nc3,along=1)  

##write means to netCDF
nc<-open.ncdf(eofnc,write=T)
put.var.ncdf(nc,varmean,ncmean)

close.ncdf(nc)


###3. ########### Calculate EOFs
##get mean dataset
nc<-open.ncdf(eofnc,write=T)
ncmean<-get.var.ncdf(nc,varmean)
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

#add EOF to netCDF file
nc<-open.ncdf(eofnc,write=T)
dim(eof)<-c(lon,lat,yy)
put.var.ncdf(nc,vareof,eof[,,1:10],count=c(-1,-1,10))

close.ncdf(nc)



############MAPS
##change EOF dimensions
dim(eof1)<-c(lon,lat) #1st EOF, in [lon,lat] 
dim(eof2)<-c(lon,lat) #
dim(eof3)<-c(lon,lat) #

##map settings
ccsm<-open.ncdf(eofnc) #extract lat, lon vectors
londim<-get.var.ncdf(ccsm,"lon") 
latdim<-get.var.ncdf(ccsm,"lat")
close.ncdf(ccsm)

data(world2HiresMapEnv) #map contours
colour<-colorRampPalette(c("red","white","blue"), space="rgb")
colour2<-colorRampPalette(c("white","yellow","orange","red"), space="rgb")

#EOF1
png(filename="EOF1.png", width=300,height=200,units="mm",res=100) 

filled.contour(londim,latdim, eof1*100, nlevels=100, zlim=c(-4,4), color.palette=colour, plot.title=title(main="1st EOF (CCSM4 Mid Holocene)", font=2, lwd=10), plot.axes={maps::map(database="world2Hires", interior=T, add=T, lwd=3)},mar=map.axes())

dev.off()

#PC and variance
png("PC_Var.png", width=500,height=400,units="mm",res=100) # PC1 and variance
plot.new()
par(mfrow=c(4,1)) #2 figures in 2 rows, 1 column
plot.ts(pc1,xlab="Years",ylab="",main="PC1",cex=3)
plot.ts(pc2,xlab="Years",ylab="",main="PC2",cex=3)
plot.ts(pc1,xlab="Years",ylab="",main="PC3",cex=3)
plot(variance[1:10]*100, main="Variances",xlab="EOF",ylab="Variance (%)",type="l")
dev.off()

#Correlation1
png(filename="Correlation1.png", width=300, height=200, units="mm", res=100) 
filled.contour(londim, latdim, corr1, nlevels=100, zlim=c(-1,1), color.palette=colour, plot.title=title(main="1st EOF Correlation", font=2, lwd=10),plot.axes={maps::map(database="world2Hires",interior=T,add=T, lwd=3)})
dev.off()

#Variance explained1
png(filename="Variance1.png", width=300, height=200, units="mm", res=100) ## variance explained1
filled.contour(londim, latdim, varexp1, nlevels=100, zlim=c(0,100), color.palette=colour2, plot.title=title(main="Variance Explained by EOF1"),plot.axes={maps::map(database="world2Hires",interior=T,add=T)})
dev.off()



##Borneo maps
#EOF
png(filename="IndoEOF1.png", width=300,height=200,units="mm",res=100) 

filled.contour(londim,latdim, eof1*100, nlevels=100, ylim=c(-20,20), xlim=c(80,150),zlim=c(-4,4), color.palette=colour, plot.title=title(main="1st EOF (CCSM4 Historical)", font=2, lwd=10), plot.axes={maps::map(database="world2Hires", interior=T, add=T, lwd=3)},mar=map.axes())

dev.off()



















