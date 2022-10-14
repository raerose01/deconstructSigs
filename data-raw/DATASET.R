## code to prepare `DATASET` dataset goes here

load("data-raw/signatures.nature2013.rda")
load("data-raw/signatures.cosmic.rda")
load("data-raw/signatures.dbs.cosmic.v3.may2019.rda")
load("data-raw/signatures.exome.cosmic.v3.may2019.rda")
load("data-raw/signatures.genome.cosmic.v3.may2019.rda")
usethis::use_data(
    signatures.nature2013,
    signatures.cosmic,
    signatures.dbs.cosmic.v3.may2019,
    signatures.exome.cosmic.v3.may2019,
    signatures.genome.cosmic.v3.may2019,
    overwrite = TRUE
)
ncol(signatures.dbs.cosmic.v3.may2019) # 78L
ncol(signatures.exome.cosmic.v3.may2019) # 96L
ncol(signatures.genome.cosmic.v3.may2019) # 96L

load("data-raw/dbs_possible.rda")
load("data-raw/tri.counts.genome.rda")
load("data-raw/tri.counts.exome.rda")
tri.counts.genome <- structure(
    tri.counts.genome[[1L]],
    names = rownames(tri.counts.genome)
)
tri.counts.exome <- structure(
    tri.counts.exome[[1L]],
    names = rownames(tri.counts.exome)
)

usethis::use_data(
    dbs_possible,
    tri.counts.genome,
    tri.counts.exome,
    internal = TRUE,
    overwrite = TRUE
)

if (!dir.exists("inst/extdata")) {
    dir.create("inst/extdata", recursive = TRUE)
}
load("data-raw/randomly.generated.tumors.rda")
saveRDS(
    randomly.generated.tumors,
    paste0("inst/extdata/", "randomly.generated.tumors", ".rds"),
    compress = "gzip"
)
load("data-raw/sample.mut.ref.rda")
saveRDS(
    sample.mut.ref,
    paste0("inst/extdata/", "sample.mut.ref", ".rds"),
    compress = "gzip"
)
load("data-raw/example.output.rda")
saveRDS(
    example.output,
    paste0("inst/extdata/", "example.output", ".rds"),
    compress = "gzip"
)
