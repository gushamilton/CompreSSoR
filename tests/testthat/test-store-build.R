test_that("GRCh37 and GRCh38 identity tables have explicit offsets", {
  for (build in c("GRCh37", "GRCh38")) {
    lengths <- CompreSSoR:::compressor_chromosome_lengths(build)
    offsets <- CompreSSoR:::compressor_chromosome_offsets(build)
    expect_identical(names(lengths), c(as.character(1:22), "X", "Y"))
    expect_identical(names(offsets), names(lengths))
    expect_equal(offsets[[1L]], 0)
    expect_equal(unname(offsets[-1L]), head(cumsum(as.numeric(lengths)), -1L))
    expect_true(all(diff(offsets) > 0))
  }
  expect_false(identical(
    CompreSSoR:::compressor_chromosome_lengths("GRCh37"),
    CompreSSoR:::compressor_chromosome_lengths("GRCh38")
  ))
})

test_that("identity encode/decode round-trips GRCh37 and GRCh38 boundaries", {
  for (build in c("GRCh37", "GRCh38")) {
    lengths <- CompreSSoR:::compressor_chromosome_lengths(build)
    input <- data.frame(
      chromosome = c("1", "chr1", "2", "chr23", "24"),
      position = c(1, lengths[["1"]], 1, lengths[["X"]], lengths[["Y"]]),
      ref = c("A", "C", "G", "T", "A"),
      alt = c("C", "T", "A", "A", "G"),
      stringsAsFactors = FALSE
    )
    encoded <- CompreSSoR:::compressor_encode_variant_identity(
      input$chromosome, input$position, input$ref, input$alt, build = build
    )
    decoded <- CompreSSoR:::compressor_decode_variant_identity(
      encoded$global_position, encoded$substitution, build = build
    )
    expect_identical(decoded$build, build)
    expect_equal(decoded$chromosome, c("1", "1", "2", "X", "Y"))
    expect_equal(decoded$position, input$position)
    expect_equal(decoded$reference_allele, input$ref)
    expect_equal(decoded$alternate_allele, input$alt)
    expect_equal(decoded$global_position, encoded$global_position)
    expect_equal(decoded$substitution, encoded$substitution)
  }
})

test_that("identity validation rejects unsupported rows and out-of-range positions", {
  lengths <- CompreSSoR:::compressor_chromosome_lengths("GRCh37")
  expect_error(
    CompreSSoR:::compressor_encode_variant_identity(
      "1", lengths[["1"]] + 1, "A", "G", build = "GRCh37"
    ),
    "outside its primary chromosome"
  )
  expect_error(
    CompreSSoR:::compressor_encode_variant_identity(
      "MT", 1, "A", "G", build = "GRCh37"
    ),
    "chromosome must be one of"
  )
  expect_error(
    CompreSSoR:::compressor_encode_variant_identity(
      "1", 1, "AT", "G", build = "GRCh37"
    ),
    "distinct single A/C/G/T"
  )
  expect_error(
    CompreSSoR:::compressor_decode_variant_identity(0, 0, build = "GRCh37"),
    "REF=ALT"
  )
  expect_error(
    CompreSSoR:::compressor_decode_variant_identity(
      sum(as.numeric(lengths)), 1, build = "GRCh37"
    ),
    "outside the selected build"
  )
})

test_that("identity manifest records build-specific provenance", {
  manifest <- CompreSSoR:::compressor_identity_manifest(
    input_build = "hg19", stored_build = "GRCh37",
    reference = list(id = "reference-37", sha256 = "ref-hash"),
    chain = list(id = "chain-37-38", sha256 = "chain-hash")
  )
  expect_equal(manifest$schema, "compressor_variant_identity_v1")
  expect_equal(manifest$input_build, "GRCh37")
  expect_equal(manifest$stored_build, "GRCh37")
  expect_equal(manifest$assembly_identifier, "GRCh37")
  expect_equal(manifest$chromosome_table$id, "grch37_primary_1_22_X_Y")
  expect_equal(manifest$chromosome_table$version, 1L)
  expect_equal(manifest$reference$id, "reference-37")
  expect_equal(manifest$chain$direction, NULL)
  expect_equal(manifest$chain$sha256, "chain-hash")

  chain_path <- tempfile(fileext = ".chain.gz")
  writeLines("chain fixture", chain_path)
  on.exit(unlink(chain_path), add = TRUE)
  hashed <- CompreSSoR:::compressor_identity_manifest(
    input_build = "GRCh37", stored_build = "GRCh38", chain = chain_path
  )
  expect_equal(hashed$chain$direction, "GRCh37-to-GRCh38")
  expect_equal(hashed$chain$sha256,
               digest::digest(chain_path, algo = "sha256", file = TRUE))
  expect_false("path" %in% names(hashed$chain))

  cross_build <- CompreSSoR:::compressor_identity_manifest(
    input_build = "GRCh37", stored_build = "GRCh38",
    chain = list(id = "chain-37-38", sha256 = "chain-hash")
  )
  expect_equal(cross_build$chain$direction, "GRCh37-to-GRCh38")
  expect_true(cross_build$chain$required)
})
