# 1. Get the list of packages you have installed, 
#    use priority to exclude base and recommended packages.
#    that may have been distributed with R.
pkgList <- installed.packages(priority='NA')[,'Package']

# 2. Find out which packages are on CRAN and R-Forge.  Because
#    of R-Forge build capacity is currently limiting the number of
#    binaries available, it is queried for source packages only.
CRANpkgs <- available.packages(
  contriburl=contrib.url('http://cran.r-project.org'))[,'Package']
forgePkgs <- available.packages(
  contriburl=contrib.url('http://R-forge.R-project.org', type='source')
)[,'Package']

# 3. Calculate the set of packages which are installed on your machine,
#    not on CRAN but also present on R-Force.
pkgsToUp <- intersect(setdiff(pkgList, CRANpkgs), forgePkgs)

# 4. Update the packages, using oldPkgs to restrict the list considered.
update.packages(checkBuilt=TRUE, ask=FALSE,
  repos="http://R-forge.R-project.org",
  oldPkgs=pkgsToUp)
