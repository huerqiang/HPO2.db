\name{HPOmetadata}
\alias{HPOmetadata}
\title{Annotation of HPO Identifiers to HPO Terms}
\description{
    This data set gives mappings between HPO identifiers and their
    respective terms.
}
\details{
    Each HPO identifier is mapped to a \code{HPOTerms} object that has 2 slots:
    name: HPO name;
    value: The value


    Mappings were based on data provided by: Human Phenotype Ontology
    With a date stamp from the source of: 20230405


}

\examples{
    # Convert the object to a data.frame
    library(AnnotationDbi)
    xx <- toTable(HPOmetadata)
}
\references{
\url{https://hpo.jax.org/app/data/ontology}
}
\keyword{datasets}

