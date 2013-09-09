# Lempel-Ziv coding of a String of char "Str" with an initial alphabet "A"

lempel_ziv_entropy = function(A, Str) {
  
  
  L   = length(A)             # Length of alphabet
  A   = as.character(A)
  Str = as.character(Str)     # convert to characters
  codebook = A                # Array for codebook
  codeDec = c()
  codeBin = c()
  N       = length(Str)
  
  # some indexes
  i = 1;
  k = 1;
  
  while (i <= length(Str)) {        # For each char
    flag = 0;                       # costraint
    search = Str[i];                # symbol for research
    while (flag==0) {
      index = which(codebook == search)
      if (length(index) > 0) {               # Symbol Found
        codeDec[k] = index                   # Add to code
        i = i + 1          
        if (i <= length(Str)) {   # If the string isn't finished
          search = paste(search, Str[i], sep='');
        } else flag=1;            # exit from while
      } else {                    # Symbol ~Found
        flag = 1;                 # exit from while
        codebook[length(codebook)+1] = search;   # add it to codebook
      }
    }        
    k = k + 1;    # number of symbols of code
  }
  codeDec[is.na(codeDec)] = 0
  codeBin = paste(sapply(codeDec, function(x) dec2bin(x,ceiling(log2(max(codeDec))))), collapse='')
  entropy = -log2( nchar(codeBin) / 64*(length(Str)))
  return(list(dec = codeDec, bin = codeBin, entropy = entropy))
}

# converts integer x to binary with at least n bits
dec2bin = function(x, n){
  y <- intToBits(x)
  z <- paste((0:1)[1+ (rev(y) == 1)], collapse="")
  z1 = gsub("^0*", "", z)   
  if (length(z1) < n) z2 = paste(paste(rep("0",n-length(z1)), collapse=''), z1, sep='') else z2 = z1
  return(z2)
}

# function to compute random entropy
rand_entropy = function(sequence) {
  return(log2(length(unique(sequence))))
}

# function to compute uncorrelated entropy (Shannon)
uncorr_entropy = function(sequence) {
  tab = table(sequence)
  p   = tab / sum(tab)
  return( -sum(p * log2(p)) )
}

gzip_entropy = function(sequence, file = 'tmp.dat'){
  # write to temporary file
  to.write = file(file, "wb")
  writeBin(sequence, to.write)
  close(to.write)
  initial = file.info(file)$size
  # compute entropy
  system(paste("gzip",file))
  compres = file.info(paste(file, "gz", sep='.'))$size
  system(paste("rm", paste(file, "gz", sep='.')))
  entropy = compres / initial
  return(entropy)
}

# function to compute maximum predictability for a given entropy value
compute_predictability = function(S, N){
  predict.fano = function(PI) {
    if (PI == 0) return(-S + log2(N))
    if (PI == 1) return(-S)
    ret = - S - (PI * log2(PI) + (1 - PI) * log2(1 - PI)) + (1-PI) * log2(N-1)
    return(ret)
  }
  PI = try(uniroot(f = predict.fano, interval = c(0,1)))
  if (class(PI) == 'try-error') res = NA else res = PI$root
  return(res)
}

