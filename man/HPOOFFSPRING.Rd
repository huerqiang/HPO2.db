\name{HPOOFFSPRING}
\alias{HPOOFFSPRING}
\title{Annotation of HPO Identifiers to their Offspring} 
\description{
    This data set describes associations between HPO 
    terms and their offspring    terms, based on the directed acyclic
    graph (DAG) defined by the Human Phenotype Ontology Consortium. The format is an R
    object mapping the HPO    terms to all offspring terms, where an
    ancestor term is a more specific HPO term that is preceded
    by the given HPO term in the DAG (in other words, the children and all
    their children, etc.).
}
\details{
    Each HPO term is mapped to a vector of offspring HPO    terms.
    Mappings were based on data provided by: Human Phenotype Ontology
    With a date stamp from the source of: 20230405


}

\examples{
    # Convert the object to a list
    xx <- as.list(HPOOFFSPRING)
}
\references{
\url{https://hpo.jax.org/app/data/ontology} 
}
\keyword{datasets}

