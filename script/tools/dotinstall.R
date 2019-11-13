require(tools)

##' dotinstall
##' @rdname dotinstall
##' @name dotinstall
##' 
NULL


as.cvec <- function(x) {
  if(length(x)>1) return(x)
  if(!any(grepl(",",x))) return(x)
  x <- gsub("\\s*", "",x)
  unlist(strsplit(x,","))
}

##' Ignore lib directory.
##' 
##' @param root project root directory
##' @param lib installed packages directory
##' @param ci run svn ci?
##' @rdname ignore_libs
##' @export
.ignore_libs <- function(root=getwd(),lib="lib", ci=FALSE) {
  
  if(!missing(root) & file.exists(root)) {
    lib <- file.path(root,"lib")
  }
  if(!file.exists(lib)) stop("Could not find lib directory")
  libs <- list.files(lib, full.names=FALSE)
  libs <- c(libs, "ignore.txt", "PACKAGES", "PACKAGES.gz")
  writeLines(con=file.path(lib,"ignore.txt"), libs)
  setwd(lib)
  system("svn propset svn:ignore -F ignore.txt .")
  setwd("..")
  if(ci) system("svn ci -m \"ignoring libs\" .")
}


.write_packages <- function(pkg,...) write_PACKAGES(pkg,...)


##' Get available packages.
##' 
##' @param pkg zipped packages directory
##' @param which column names to return from \code{\link{available.packages}}
##' @param write.first call \code{\link{write_PACKAGES}} first?
##' @param ... passed to \code{\link{write_PACKAGES}}
##' @rdname get_available
##' @export 
.get_available <- function(pkg,which=c("Package"),write.first=TRUE,...) {
  pkg <- normalizePath(pkg)
  if(write.first) .write_packages(pkg,...)
  f <- available.packages(contriburl=.my_cran(pkg))[,which]
  return(f)
}
.my_cran <- function(pkg) {
  paste0("file://", pkg) 
}


##' Install locally archived packages.
##' 
##' @param ... quoted names of local packages to install
##' @param root project root directory
##' @param lib installed packages directory
##' @param pkg zipped packages directory
##' @param verbose not used
##' @rdname install_local
##' @export
.install_local <- function(...,root=getwd(),lib="lib", pkg="pkg",verbose=FALSE,type="source") {
  
  pkgs <- list(...)
  pkgs <- unlist(lapply(pkgs, as.cvec))
  
  scriptDir <- root
  pkg <- normalizePath(file.path(scriptDir, pkg))
  lib <- normalizePath(file.path(scriptDir, lib))
  
  if(!file.exists(pkg)) stop("Could not find pkg directory")
  if(!file.exists(lib)) stop("Could not find lib directory")
  
  home <- Sys.getenv("HOME")
  Sys.setenv(HOME="")
  on.exit(Sys.setenv(HOME=home))
  
  .libPaths(lib)
  require(tools)
  
  f <- .get_available(pkg)
  
  N <- length(f)
  
  if(N==0) stop(paste("Could not find any packages in",pkg))
  
  cat("\n\nFound ", N, " packages in local repository\n")
  
  if(is.character(pkgs)) {
    f <- f[f %in% pkgs]
    if(length(f)==0) stop("Could not find requested packages: ", pkgs, "...")
    cat("Installing ", paste(pkgs, collapse=" "), "\n\n")
  } else {
    cat("Installing all packages in local repository...\n\n")
  }
  install.packages(f,contriburl=.my_cran(pkg),type=type, lib="lib", INSTALL_opts="--no-multiarch")
}


##' Install packages from CRAN or other CRAN mirror.
##' 
##' @param ... character names of R packages to install
##' @param root project root directory
##' @param repos repository from which to get packages
##' @param lib installed packages directory
##' @param pkg zipped packages directory
##' @param pkg.win directory for storing win.binary versions of packages
##' @param get.win logical; if \code{TRUE}, win.binary packages will be retreived
##' @rdname install_cran
##' @export
##' @details
##' The R-Forge repository location \code{http://R-Forge.R-project.org} is automatically appended to \code{repos}.
.install_cran <- function(..., root=getwd(),repos="http://cran.us.r-project.org",dep=TRUE,
                          lib="lib", pkg="pkg",type="source", pkg.win=file.path(pkg,"win"), get.win=FALSE) {
  
  pkgs <- list(...)
  pkgs <- unlist(lapply(pkgs, as.cvec))
  
  repos <- c(repos, "http://R-Forge.R-project.org")
  
  if(get.win & !file.exists(pkg.win)) stop(paste("Windows binary directory", pkg.win, "doesn't exist."))
  
  
  pkg <- normalizePath(path.expand(file.path(root,pkg)))
  lib <- normalizePath(path.expand(file.path(root,lib)))
  if(!file.exists(pkg)) stop("Could not find pkg directory ", pkg)
  if(!file.exists(lib)) stop("Could not find lib directory ", lib)
  
  cat("Installing: ", paste( pkgs, collapse=", "), "\n\n")
  cat("library location: ", lib, "\n")
  cat("destdir: ", pkg, "\n\n")
  
  home <- Sys.getenv("HOME")
  Sys.setenv(HOME="")
  on.exit(Sys.setenv(HOME=home))
  
  if(length(pkgs) <1) stop("Please specify packages to install from CRAN")
  
  .libPaths(lib)
  require(tools)
  
  pkgs <- .pkg_list(pkgs)
  
  if(length(pkgs$dep)>0 & dep) {
    cat("dotinstall found these dependencies:\n")
    cat(paste("   ", pkgs$dep), sep="\n")
    cat("\n")
  }
  if(!dep) pkgs$dep <- character(0)
  pkgs <- sort(unique(unlist(pkgs)))  
  
  install.packages(pkgs,repos=repos, type=type, destdir=pkg,lib=lib,INSTALL_opts="--no-multiarch")
  
  if(get.win) {
    cat("\n\nDownloading windows binaries ...")
    if(!file.exists(pkg.win)) stop(paste("Directory ", pkg.win, " doesn't exist"))
    .get_win_binaries(pkgs, destdir=pkg.win, repos=repos)
  }
  return(invisible(NULL))
}

##' Return a matrix of locally-available packages.
##' 
##' @param root root project directory 
##' @param pkg zipped packages directory
##' @param ... passed along
##' @rdname local_repos
##' @export
.local_repos <- function(root=getwd(), pkg="pkg",...) {
  
  if(!missing(root) & file.exists(root)) {
    pkg <- file.path(root,pkg)
  } else {
    scriptDir <- root
    pkg <- file.path(scriptDir, pkg)
  }
  
  if(!file.exists(pkg)) stop("Could not find pkg directory")
  
  require(tools)
  
  f <- .get_available(pkg, which=c("Package", "Version"))
  
  N <- nrow(f)
  
  cat("\nFound ", N, " packages in local repository")
  return(f)
  
}


##' Get a matrix listing locally-installed packages.
##' 
##' @param root project root directory
##' @param lib installed packages directory
##' @param ... passed along
##' @rdname local_installed
##' @export
.local_installed <- function(root=getwd(),lib="lib",...) {
  if(!missing(root) & file.exists(root)) {
    lib <- file.path(root,lib)
  } else {
    scriptDir <- root
    lib <- file.path(scriptDir, lib)
  }
  
  if(!file.exists(lib)) stop("Could not find lib directory")
  
  .libPaths(lib)
  library(tools)
  f <- installed.packages(lib)[,c("Package", "Version")]
  
  N <- nrow(f)
  
  cat("\nFound ", N, " packages installed locally\n\n")
  print(f)
  return(invisible(as.data.frame(f,stringsAsFactors=FALSE)))
  
}

.pkg_list <- function(x,base=FALSE,...) {
  force(x)
  x <- .no_base(x,...)
  y <- .get_dependencies(x,...)
  if(!base) {
    x <- sort(unique(.no_base(x,...)))
    y <- sort(unique(.no_base(y,...)))
  }
  return(list(pkg=x, dep=y))
}
##' Look up package dependencies.
##' 
##' @param x character vector of package names
##' @param which passed to \code{\link{package_dependencies}}
##' @param db passed to \code{\link{package_dependencies}}
##' @rdname get_dependencies
##' @export
.get_dependencies <- function(x, which=c("Depends", "Imports", "LinkingTo"), db=NULL,...) {
  require(tools)
  if(is.null(db)) db <- utils::available.packages()
  x <- .no_base(x)
  return(tools::package_dependencies(x, which=which,recursive=TRUE, db=db))
}

.get_win_binaries <- function(x,destdir='.',base=FALSE,...) {
  if(base) x <- .no_base(x)
  utils::download.packages(x, destdir=destdir, type="win.binary",...) 
}

.no_base <- function(x,...) {
  y <- names(which(utils::installed.packages()[, "Priority"] == "base"))
  setdiff(unlist(x), y) 
}

.get_audited <- function(pkg,lib,install=TRUE) {
  dest <- normalizePath(file.path(pkg,"audited_1.9.tar.gz"))
  download.file("http://metrumrg-soft.s3.amazonaws.com/audited/audited_1.9.tar.gz",
                destfile=dest) 
  if(install) install.packages(dest, lib=lib,type="source", repos=NULL, INSTALL_opts="--no-multiarch")
}


