\name{HPOANCESTOR}
\alias{HPOANCESTOR}
\title{Annotation of HPO Identifiers to their Ancestors} 
\description{
    This data set describes associations between HPO 
    terms and their ancestor    terms, based on the directed acyclic
    graph (DAG) defined by the Human Phenotype Ontology Consortium. The format is an R
    object mapping the HPO    terms to all ancestor terms, where an
    ancestor term is a more general HPO term that precedes
    the given HPO term in the DAG (in other words, the parents, and all
    their parents, etc.).
}
\details{
    Each HPO    term is mapped to a vector of ancestor HPO    terms.
    Mappings were based on data provided by: Human Phenotype Ontology
    With a date stamp from the source of: 20230405


}

\examples{
    # Convert the object to a list
    xx <- as.list(HPOANCESTOR)
}

\references{
\url{https://hpo.jax.org/app/data/ontology} 
}
\keyword{datasets}

