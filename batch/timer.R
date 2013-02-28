tic <- function(name='default', gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self"))
{
  type <- match.arg(type)
  assign(".type", type, envir=baseenv())
  if(gcFirst) gc(FALSE)
  tic = list()
  tryCatch({tic <- get(".tic", envir=baseenv())},
           error=function(e){  })
  tic[name] <- proc.time()[type]         
  assign(".tic", tic, envir=baseenv())
  invisible(tic)
}

toc <- function(name='default',prefixStr=NA)
{
  type <- get(".type", envir=baseenv())
  tic  <- get(".tic", envir=baseenv())
  dt   <- proc.time()[type] - as.numeric(tic[name])
  
  # must be ints...
  d <- floor(dt / 3600 / 24)
  h <- floor(dt / 3600) %% 24
  m <- floor(dt / 60)   %% 60
  s <- dt %% 60
  #f <- s - floor(s)
  #s <- floor(s)
  if(is.na(prefixStr)) prefixStr <- name
  print(paste(prefixStr,': ',sprintf('%02i:%02i:%02i:%05.2f',d,h,m,s),sep=''))

}
