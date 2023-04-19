\name{HPOCHILDREN}
\alias{HPOCHILDREN}
\title{Annotation of HPO Identifiers to their Children} 
\description{
This data set describes associations between HPO 
    terms and their direct children    terms, based on the directed acyclic
    graph (DAG) defined by the Human Phenotype Ontology Consortium. The format is an R
    object mapping the HPO    terms to all direct children terms, where a
    direct child term is a more specific HPO term that is immediately preceded
    by the given HPO term in the DAG.
}
\details{
    Each HPO term is mapped to a vector of children HPO terms.
    Mappings were based on data provided by: Human Phenotype Ontology
    With a date stamp from the source of: 20230405


}
\examples{
    # Convert the object to a list
    xx <- as.list(HPOCHILDREN)
}


\references{
\url{https://hpo.jax.org/app/data/ontology} 
}
\keyword{datasets}

