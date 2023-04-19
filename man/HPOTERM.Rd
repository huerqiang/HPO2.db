\name{HPOTERM}
\alias{HPOTERM}
\title{Annotation of HPO Identifiers to HPO Terms}
\description{
    This data set gives mappings between HPO identifiers and their
    respective terms.
}
\details{
    Each HPO identifier is mapped to a \code{HPOTerms} object that has 2 slots:
    do_id: HPO Identifier;
    Term: The term for that HPO id


    Mappings were based on data provided by: Human Phenotype Ontology
    With a date stamp from the source of: 220230405


}

\examples{
    # Convert the object to a list
    xx <- as.list(HPOTERM)
}
\references{
\url{https://hpo.jax.org/app/data/ontology}
}
\keyword{datasets}

