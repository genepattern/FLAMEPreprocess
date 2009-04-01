unzip.file <- function(input.filename, temp.dir){
   isWindows <- Sys.info()[["sysname"]]=="Windows"

   if(isWindows)
   {
       zip.unpack(input.filename, dest = temp.dir)
   }
   else
   {
       zip <- getOption("unzip")
       system(paste(zip, " -j -q '", input.filename, "' -d ", temp.dir, sep=''))
   }
}