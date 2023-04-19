\name{HPOPARENTS}
\alias{HPOPARENTS}
\title{Annotation of HPO Identifiers to their Parents} 
\description{
    This data set describes associations between HPO
    terms and their direct parent    terms, based on the directed acyclic
    graph (DAG) defined by the Human Phenotype Ontology Consortium. The format is an R
    object mapping the HPO    terms to all direct parent terms, where a
    direct parent term is a more general HPO term that immediately precedes
    the given HPO term in the DAG.
}
\details{
    Each HPO    term is mapped to a named vector of HPO    terms. The name
    associated with the parent term will be either \emph{isa}, \emph{partof},
    where \emph{isa} indicates that the child term is a more specific version
    of the parent, and    \emph{partof} indicate that the
    child term is a part of the parent. 


    Mappings were based on data provided by: Human Phenotype Ontology
    With a date stamp from the source of: 220230405


}
\examples{
    xx <- as.list(HPOPARENTS)
}
\references{
\url{https://hpo.jax.org/app/data/ontology} 
}
\keyword{datasets}

