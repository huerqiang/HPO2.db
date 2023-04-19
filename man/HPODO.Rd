\name{HPODO}
\alias{HPODO}
\title{Annotation of  HPO Identifiers to DO ID} 
\description{
    This data set describes associations between HPO 
    ids and do ids, based on the Human Phenotype Ontology. 
    The format is an R
    object mapping the HPO ids to do ids.
}
\details{
    Each HPO id is mapped to a vector of do ids.
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
    xx <- as.list(HPODO)
}

\references{
\url{https://github.com/mapping-commons/mh_mapping_initiative}
}
\keyword{datasets}
