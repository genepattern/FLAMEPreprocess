#LoadPackages <- function(package.names)
#{
#    #source("http://bioconductor.org/biocLite.R")
#    for (i in 1:length(package.names))
#    {
#	    package.name = package.names[i]
#	    installed <- installed.packages()[,1]
#	    if (!any(installed == package.name))
#	    {
#		    #install.packages(package.name, repos = "http://cran.r-project.org")
#		    biocLite(package.name)
#	    }
#    }
#}

cleanup <- function()
{
    setwd(wkdir)

    temp.files <- list.files(temp.dir, full.names=TRUE)
    unlink(temp.files, recursive = TRUE)     
    unlink(temp.dir, recursive = TRUE)

    files <- list.files(all.files=TRUE)
    for (i in 1:length(files))
    {
        if(regexpr(paste(".zip","$",sep=""), tolower(files[[i]]))[[1]] == -1
            && tolower(files[[i]]) != "stderr.txt" && tolower(files[[i]]) != "cmd.out"
            && tolower(files[[i]]) != "stdout.txt")
        {
            file.remove(files[[i]])
        }
    }
}

parseCmdLine <- function(...)
{
    suppressMessages(preprocess(...))
}

preprocess <- function(
libdir, #full pathname to where FLAME program files are kept
dataset, #full pathname to zipped data files
remove_dead, #whether to remove dead cells automatically
channels,
channel_names,
scatter_channels = "1,2",
filetype, #{fcs, txt}
transformation = "logicle", #none, logicle, arcsinh
r, #10000,262144
logicle_cofactor = 3,
arcsinh_cofactor = 250,
output_prefix #<studyname_dataname>
){

zip.ext <- regexpr(paste(".zip","$",sep=""), tolower(dataset))
if(zip.ext[[1]] == -1)
{
    stop("Input file must be of type zip ")
}

#gather channel info
columns <- as.numeric(strsplit(channels,',')[[1]])
if(length(na.omit(columns)) != length(columns))
{
    stop("channels must be integer values ")
}
print(columns)

scatter_channels <- as.numeric(strsplit(scatter_channels,',')[[1]])

headers <- strsplit(channel_names,',')[[1]]
print(headers)
dim = length(columns)

if(length(columns) != length(headers))
{
    stop(paste("The number of channels:", length(columns),
        "is not equal to the number of channel names:", length(headers), " "))
}

source(paste(libdir,"common.R",sep='/'))
source(paste(libdir,"zip.R",sep='/'))
source(paste(libdir,"unzip.R",sep='/'))
source(paste(libdir,"EmSkew.R",sep='/'))
source(paste(libdir,"transformfcs.R",sep='/'))

if(libdir!='')
{
    setLibPath(libdir)
    install.required.packages(libdir)
}

suppressMessages(library(flowCore))

lTrans <<- logicleTransform("logicle", d=logicle_cofactor,r=10000)
lTrans.2 <<- logicleTransform("logicle", d=logicle_cofactor, r=262144)
arcsinh <- function(x,c=arcsinh_cofactor) {
	f = log(x/c + sqrt((x/c)^2+1))
	return(f)
}

parDefault <- function(exp){
  vm <- data.frame(labelDescription=c(name="Name of Parameter",
                   desc="Description of Parameter",
                   range="Range of Parameter",
                   minRange="Minimum Parameter Value after Transformation",
                   maxRange="Maximum Parameter Value after Transformation"))
  pd <- data.frame(name=colnames(exp), desc=colnames(exp),
                   range=apply(exp, 2, max, na.rm=TRUE),
                   minRange=apply(exp, 2, min, na.rm=TRUE),
                   maxRange=apply(exp, 2, max, na.rm=TRUE))
  new("AnnotatedDataFrame", pd, vm)
}

wkdir <<- getwd()

setwd(libdir)

isWindows <- Sys.info()[["sysname"]]=="Windows"
if(isWindows)
{
    file.copy("emskew.dll", to = wkdir)
}
else
{
    file.copy("emskew.so", to = wkdir)
}
setwd(wkdir)

temp.dir <<- paste(wkdir, "temp", sep="/")
dir.create(temp.dir)
on.exit(cleanup())

unzip.file(dataset, temp.dir)
datafiles <- list.files(temp.dir, full.names = TRUE)

#filter, log
if (filetype == "fcs") {
	filenames<-c()
	for (i in 1:length(datafiles)) {
		cat(i,'\n')
		filename <- strsplit(datafiles[i], "\\.fcs")
		file<-read.FCS(datafiles[i])
		data <- exprs(file)

		if (remove_dead == "T") {
#			remove dead cells using skew-t g=4
			scatters <- data[,1:2]
			obj <- EmSkew(dat=scatters,ng=4,dist=4)
			FSCmus <- obj$mu[1,]
			SSCmus <- obj$mu[2,]
			dead <- which(FSCmus == min(FSCmus))
			clusters <- obj$clust
			dead.cells <- scatters[which(clusters == dead),]

			if (.Platform$OS.type == "windows")
            {
			    png(paste(filename,"remove_dead.png",sep='.'),height=960,width=960)
			}
			else
			{
                library(Cairo)
			    CairoPNG(paste(filename,"remove_dead.png",sep='.'),height=960,width=960)
			}

			plot(scatters,pch='.',main = paste(filename,"remove_dead",sep='.'))
			points(dead.cells,pch='.',col='red')
			dev.off()

			data <- data[which(clusters != dead),]
			#dead-removed data
		}

		logicledata <- c()
		for (c in columns) {
			thiscolumn <- matrix(data[,c],dimnames = list(c(1:(dim(data)[1])),c("this")))
			if (transformation=="logicle") {
				if (!any(c==scatter_channels)) {#if (c!=1 && c!=2) {
					param <- parDefault(thiscolumn)
					ff <- new("flowFrame",exprs=as.matrix(thiscolumn),parameters=param)
					if (r==10000) {
						transform <- transform(ff, `this` = lTrans(as.real(`this`)))
					}
					if (r==262144) {
						transform <- transform(ff, `this` = lTrans.2(as.real(`this`)))
					}
					thiscolumn <- exprs(transform)
				}
			}
			if (transformation=="arcsinh") {
				if (!any(c==scatter_channels)) {#if (c!=1 && c!=2) {
					thiscolumn <- arcinh(thiscolumn)
				}
			}
			logicledata <- cbind(logicledata,thiscolumn)
		}
		logicledata <- na.omit(logicledata)
		write.table(logicledata, file = paste(filename,".preprocessed.txt", sep = ""),sep = "\t", row.names = F,  col.names = headers, quote = F)
	}
}

if (filetype == "txt") {
	for (i in 1:length(datafiles)) {
		cat(i,'\n')
		filename <- strsplit(datafiles[i], "\\.txt")
		file <- read.table(datafiles[i],header=F,skip=1)
		data<- data.frame(file)
#		data <- subset(file, select = columns)

		if (remove_dead == "T") {
#			remove dead cells using skew-t g=4
			scatters <- data[,1:2]
			obj <- EmSkew(dat=scatters,ng=4,dist=4)
			FSCmus <- obj$mu[1,]
			SSCmus <- obj$mu[2,]
			dead <- which(FSCmus == min(FSCmus))
			clusters <- obj$clust
			dead.cells <- scatters[which(clusters == dead),]

			if (.Platform$OS.type == "windows")
            {
			    png(paste(filename,"remove_dead.png",sep='.'),height=960,width=960)
			}
			else
			{
                library(Cairo)
			    CairoPNG(paste(filename,"remove_dead.png",sep='.'),height=960,width=960)
			}

			plot(scatters,pch='.',main = paste(filename,"remove_dead",sep='.'))
			points(dead.cells,pch='.',col='red')
			dev.off()

			data <- data[which(clusters != dead),]
			#dead-removed data
		}

		logicledata <- c()#matrix(nrow= nrow(data),ncol=length(columns))
		for (c in columns) {
			thiscolumn <- matrix(data[,c],dimnames = list(c(1:(dim(data)[1])),c("this")))
			if (transformation=="logicle") {
				if (!any(c==scatter_channels)) {#if (c!=1 && c!=2) {
					param <- parDefault(thiscolumn)
					ff <- new("flowFrame",exprs=as.matrix(thiscolumn),parameters=param)
					if (r==10000) {
						transform <- transform(ff, `this` = lTrans(as.real(`this`)))
					}
					if (r==262144) {
						transform <- transform(ff, `this` = lTrans.2(as.real(`this`)))
					}
					thiscolumn <- exprs(transform)
				}
			}
			if (transformation=="arcsinh") {
				if (!any(c==scatter_channels)) {#if (c!=1 && c!=2) {
					thiscolumn <- arcsinh(thiscolumn)
				}
			}
			logicledata <- cbind(logicledata,thiscolumn)
		}
		logicledata <- na.omit(logicledata)
		write.table(logicledata, file = paste(filename,".preprocessed.txt", sep = ""),sep = "\t", row.names = F,  col.names = headers, quote = F)
	}
}

##ouput files as zip

setwd(temp.dir)

if(remove_dead == "T")
{
    zip.file(libdir = libdir, files = "*.remove_dead.png", outfile = paste(wkdir, "/", output_prefix,".Remove-dead.scatterplots.zip",sep='')) #change libdir to where zip.exe is saved
}

zip.file(libdir = libdir, files = "*.preprocessed.txt", outfile = paste(wkdir, "/", output_prefix,".PreprocessedData.zip",sep='')) #change libdir to where zip.exe is saved

setwd(wkdir)
}

install.required.packages <- function(libdir)
{
	if(!is.package.installed(libdir, "robustbase"))
    {
		install.package(libdir,,, "robustbase_0.4-5.tar.gz")
	}
	if(!is.package.installed(libdir, "pcaPP"))
    {
		install.package(libdir,,, "pcaPP_1.6.tar.gz")
	}
    if(!is.package.installed(libdir, "mvtnorm"))
    {
		install.package(libdir,,, "mvtnorm_0.9-4.tar.gz")
	}

	if(!is.package.installed(libdir, "rrcov"))
	{
		install.package(libdir,,, "rrcov_0.4-08.tar.gz")
	}

    #-----------------------------------------
    if(!is.package.installed(libdir, "rpanel"))
    {
		install.package(libdir,,, "rpanel_1.0-5.tar.gz")
	}
    if(!is.package.installed(libdir, "ks"))
    {
		install.package(libdir,,, "ks_1.5.10.tar.gz")
	}

	if(!is.package.installed(libdir, "feature"))
	{
		install.package(libdir,,, "feature_1.2.0.tar.gz")
	}

    #-------------------------------------------
    if(!is.package.installed(libdir, "Biobase"))
	{
		install.package(libdir,,, "Biobase_2.2.2.tar.gz")
	}
	if(!is.package.installed(libdir, "graph"))
	{
		install.package(libdir,,, "graph_1.20.0.tar.gz")
	}
    if(!is.package.installed(libdir, "flowCore"))
	{
		install.package(libdir,,, "flowCore_1.8.1.tar.gz")
	}
}