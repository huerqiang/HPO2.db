
datacache <- new.env(hash=TRUE, parent=emptyenv())

HPO <- function() showQCData("HPO", datacache)
HPO_dbconn <- function() dbconn(datacache)
HPO_dbfile <- function() dbfile(datacache)
HPO_dbschema <- function(file="", show.indices=FALSE) dbschema(datacache,
    file=file, show.indices=show.indices)
HPO_dbInfo <- function() dbInfo(datacache)

.onLoad <- function(libname, pkgname)
{
    dbfile <- system.file("extdata", "HPO.sqlite", package=pkgname,
        lib.loc=libname)
    assign("dbfile", dbfile, envir=datacache)
    dbconn <- dbFileConnect(dbfile)
    assign("dbconn", dbconn, envir=datacache)

    HPODb <- setRefClass("HPODb", contains="GODb")
    ## Create the OrgDb object
    sPkgname <- sub(".db$","",pkgname)
    txdb <- loadDb(system.file("extdata", paste(sPkgname,
        ".sqlite",sep=""), package=pkgname, lib.loc=libname),
        packageName=pkgname)
    dbObjectName <- getFromNamespace("dbObjectName", "AnnotationDbi")
    dbNewname <- dbObjectName(pkgname,"HPODb")
    ns <- asNamespace(pkgname)
    assign(dbNewname, txdb, envir=ns)
    namespaceExport(ns, dbNewname)

    ## Create the AnnObj instances
    ann_objs <- createAnnObjs.HPO_DB("HPO", "HPO", dbconn, datacache)
    mergeToNamespaceAndExport(ann_objs, pkgname)
}

.onUnload <- function(libpath)
{
    dbFileDisconnect(HPO_dbconn())
}

