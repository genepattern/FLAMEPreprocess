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
            && tolower(files[[i]]) != "stdout.txt"
            && tolower(files[[i]]) != ".epilogue.pbs"
            && tolower(files[[i]]) != "command.pbs"
            && tolower(files[[i]]) != ".command.pbs"
            && tolower(files[[i]]) != ".epilogue.sh"
)
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
transformation = "logicle", #none, logicle, arsinh or arcsinh
r, #10000=4-decade,262144=18-bit
logicle_cofactor = 3,
arsinh_cofactor = 250,
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

if(Sys.getenv("R_LIBS") != '')
{
    setLibPath(c(Sys.getenv("R_LIBS"), .libPaths()))
}

suppressMessages(library(flowCore))

logicle_cofactor <- as.numeric(logicle_cofactor)
if(is.na(logicle_cofactor))
{
    stop("logical cofactor must be a number")
}

arsinh_cofactor <- as.numeric(arsinh_cofactor)
if(is.na(arsinh_cofactor))
{
    stop("arsinh cofactor must be a number")
}


lTrans <<- logicleTransform("logicle", d=logicle_cofactor,r=10000)
lTrans.2 <<- logicleTransform("logicle", d=logicle_cofactor, r=262144)
arsinh <- function(x,c=arsinh_cofactor) {
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
datafiles.names <- list.files(temp.dir, full.names=FALSE)


#filter, log
if (filetype == "fcs") {
	filenames<-c()
	for (i in 1:length(datafiles)) {
		cat("processing file:" , datafiles.names[i],'\n')
		filename <- strsplit(datafiles[i], "\\.fcs")
		file<-read.FCS(datafiles[i])
		data <- exprs(file)
		validate.num.channels(data, columns)

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

		transformed_data <- c()
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
			if (transformation=="arsinh") {
				if (!any(c==scatter_channels)) {#if (c!=1 && c!=2) {
					thiscolumn <- arsinh(thiscolumn)
				}
			}
			transformed_data <- cbind(transformed_data,thiscolumn)
		}
		transformed_data <- na.omit(transformed_data)
		write.table(transformed_data, file = paste(filename,".preprocessed.txt", sep = ""),sep = "\t", row.names = F,  col.names = headers, quote = F)
	}
}

if (filetype == "txt") {
	for (i in 1:length(datafiles)) {
		cat("processing file:" , datafiles.names[i],'\n')
		filename <- strsplit(datafiles[i], "\\.txt")
		file <- read.table(datafiles[i],header=F,skip=1)
		data<- data.frame(file)
		validate.num.channels(data, columns)

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

		transformed_data <- c()#matrix(nrow= nrow(data),ncol=length(columns))
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
			if (transformation=="arsinh") {
				if (!any(c==scatter_channels)) {#if (c!=1 && c!=2) {
					thiscolumn <- arsinh(thiscolumn)
				}
			}
			transformed_data <- cbind(transformed_data,thiscolumn)
		}
		transformed_data <- na.omit(transformed_data)
		write.table(transformed_data, file = paste(filename,".preprocessed.txt", sep = ""),sep = "\t", row.names = F,  col.names = headers, quote = F)
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

validate.num.channels <- function(data, channels)
{
    max.channel <- max(channels)
    if(max.channel > ncol(data))
    {
        stop("Maximum channel number specified", max.channel, "is greater than number of channels in dataset", ncol(data))
    }
}

install.required.packages <- function(libdir)
{
	if(!is.package.installed(libdir, "robustbase"))
    {
		install.package(libdir, "robustbase_0.4-5.zip", "robustbase_0.4-5.tgz", "robustbase_0.4-5.tar.gz")
	}
	if(!is.package.installed(libdir, "pcaPP"))
    {
		install.package(libdir, "pcaPP_1.6.zip", "pcaPP_1.6.tgz", "pcaPP_1.6.tar.gz")
	}
    if(!is.package.installed(libdir, "mvtnorm"))
    {
		install.package(libdir, "mvtnorm_0.9-7.zip", "mvtnorm_0.9-7.tgz", "mvtnorm_0.9-7.tar.gz")
	}
	if(!is.package.installed(libdir, "rrcov"))
	{
		install.package(libdir, "rrcov_0.5-01.zip", "rrcov_0.5-01.tgz", "rrcov_0.5-01.tar.gz")
	}

    #-----------------------------------------
    if(!is.package.installed(libdir, "rpanel"))
    {
		install.package(libdir, "rpanel_1.0-5.zip", "rpanel_1.0-5.tgz", "rpanel_1.0-5.tar.gz")
	}
    if(!is.package.installed(libdir, "ks"))
    {
		install.package(libdir, "ks_1.6.5.zip", "ks_1.6.5.tgz", "ks_1.6.5.tar.gz")
	}
	if(!is.package.installed(libdir, "feature"))
	{
		install.package(libdir, "feature_1.2.3.tar.gz", "feature_1.2.3.tar.gz", "feature_1.2.3.tar.gz")
	}

    #-------------------------------------------
    if(!is.package.installed(libdir, "Biobase"))
	{
	    if(length(grep(R.version$os, "darwin9")) != 0)
	    {
	        mac_package <- "Biobase_2.2.2_leopard.tgz"
	    }
	    else
	    {
	        mac_package <- "Biobase_2.2.2_tiger.tgz"
	    }

		install.package(libdir, "Biobase_2.2.2.zip", mac_package, "Biobase_2.2.2.tar.gz")
	}
	if(!is.package.installed(libdir, "graph"))
	{
		install.package(libdir, "graph_1.22.2.zip", "graph_1.22.2.tgz", "graph_1.22.2.tar.gz")
	}
    if(!is.package.installed(libdir, "flowCore"))
	{
	    if(length(grep(R.version$os, "darwin9")) != 0)
	    {
	        mac_package <- "flowCore_1.8.3_leopard.tgz"
	    }
	    else
	    {
	        mac_package <- "flowCore_1.8.3_tiger.tgz"
	    }

		install.package(libdir, "flowCore_1.8.3.zip", mac_package, "flowCore_1.8.3.tar.gz")
	}
}