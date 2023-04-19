\name{HPOMPO}
\alias{HPOMPO}
\title{Annotation of  HPO Identifiers to MPO ID} 
\description{
    This data set describes associations between HPO 
    ids and MPO ids, based on the Human Phenotype Ontology. 
    The format is an R
    object mapping the HPO ids to MPO ids.
}
\details{
    Each HPO id is mapped to a vector of MPO ids.
    Mappings were based on data provided by: The Human Phenotype Ontology
}

\value{
    \code{HPO_dbconn}: a DBIConnection object representing an
    open connection to the package annotation DB.

    \code{HPO_dbfile}: a character string with the path to the
    package annotation DB.

    \code{HPO_dbschema}: none (invisible \code{NULL}).

    \code{HPO_dbInfo}: none (invisible \code{NULL}).
}

\examples{
    # Convert the object to a list
    xx <- as.list(HPOMPO)
}

\references{
\url{https://github.com/mapping-commons/mh_mapping_initiative}
}
\keyword{datasets}
