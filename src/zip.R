zip.file <- function(libdir, files, outfile) {
  #outfile <- check.extension(outfile, ".zip")
  isWindows <- Sys.info()[["sysname"]]=="Windows"
  if(isWindows) {
      libdir <- sub("Program Files", "Progra~1", libdir)
      libdir <- sub("Documents and Settings", "Docume~1", libdir)
      last.char <- substr(libdir, nchar(libdir), nchar(libdir))
      if(last.char != "/" && last.char!="\\") {
          libdir <- paste(libdir, "/", sep='')
      }
      system(paste(libdir, "zip.exe -q -m ", outfile, " ", files, sep=''))
  } else {
      system(paste("zip -q -m ", outfile, " ", files, sep=''))
  }
}