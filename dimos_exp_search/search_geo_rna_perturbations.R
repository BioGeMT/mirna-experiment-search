#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DBI)
  library(GEOmetadb)
  library(optparse)
  library(parallel)
  library(readr)
  library(RSQLite)
  library(stringr)
})

option_list <- list(
  make_option(
    "--targets-file",
    dest = "targets_file",
    help = "File with one mature miRNA name per line"
  ),
  make_option(
    "--mode",
    default = "both",
    help = "Perturbation mode to search: oe, ko, or both [default: %default]"
  ),
  make_option(
    "--date",
    default = "",
    help = "Optional comma-separated GEO submission years, e.g. 2023,2024"
  ),
  make_option(
    "--output",
    default = "geo_rna_perturbation_hits.tsv",
    help = "Output TSV path [default: %default]"
  ),
  make_option(
    "--ambiguous-output",
    dest = "ambiguous_output",
    default = "",
    help = "Optional TSV path for mixed_or_ambiguous rows removed from the main output"
  ),
  make_option(
    "--cpus",
    type = "integer",
    default = 1,
    help = "Number of parallel workers for --targets-file [default: %default]"
  )
)

opts <- parse_args(OptionParser(option_list = option_list))
db_path <- "GEOmetadb.sqlite"

if (is.null(opts$targets_file) || !nzchar(opts$targets_file)) {
  stop("--targets-file is required", call. = FALSE)
}
mature_mirnas <- read_lines(opts$targets_file)
mature_mirnas <- str_trim(sub("#.*$", "", mature_mirnas))
mature_mirnas <- mature_mirnas[nzchar(mature_mirnas)]
mature_mirnas <- unique(mature_mirnas)
if (length(mature_mirnas) == 0) {
  stop(
    "--targets-file must contain at least one mature miRNA name.",
    call. = FALSE
  )
}

if (!opts$mode %in% c("oe", "ko", "both")) {
  stop("--mode must be one of: oe, ko, both", call. = FALSE)
}
if (is.na(opts$cpus) || opts$cpus < 1) {
  stop("--cpus must be a positive integer", call. = FALSE)
}

date_years <- str_trim(unlist(strsplit(opts$date, ",")))
date_years <- date_years[nzchar(date_years)]
if (length(date_years) > 0 && any(!str_detect(date_years, "^\\d{4}$"))) {
  stop("--date must be a comma-separated list of four-digit years", call. = FALSE)
}
species_filter <- "Homo sapiens"

if (!file.exists(db_path)) {
  options(timeout = max(getOption("timeout"), 7200))
  downloaded <- getSQLiteFile(destdir = dirname(normalizePath(db_path, mustWork = FALSE)))
  if (normalizePath(downloaded) != normalizePath(db_path, mustWork = FALSE)) {
    file.copy(downloaded, db_path, overwrite = TRUE)
  }
}

rna_seq_terms <- c(
  "rna-seq",
  "rnaseq",
  "rna seq",
  "transcriptome sequencing",
  "expression profiling by high throughput sequencing",
  "high-throughput sequencing"
)

knockdown_terms <- c(
  "knockdown",
  "knock-down",
  "silencing",
  "sirna",
  "shrna",
  "inhibitor",
  "mirna inhibitor",
  "mir inhibitor",
  "antagomir",
  "sponge",
  "aso",
  "antisense",
  "gapmer",
  "crispri"
)

knockout_terms <- c(
  "knockout",
  "knock-out",
  "null",
  "deficient",
  "deleted",
  "deletion",
  "crispr",
  "cas9"
)

overexpression_terms <- c(
  "overexpression",
  "over-expression",
  "overexpressed",
  "over-expressed",
  "ectopic expression",
  "mimic",
  "mirna mimic",
  "mir mimic",
  "premir",
  "pre-mir"
)

selected_perturbation_terms <- switch(
  opts$mode,
  oe = overexpression_terms,
  ko = knockout_terms,
  both = c(knockout_terms, overexpression_terms)
)

text_fields <- c(
  "gse.title",
  "gse.summary",
  "gse.overall_design",
  "gse.type",
  "gsm.title",
  "gsm.source_name_ch1",
  "gsm.characteristics_ch1",
  "gsm.treatment_protocol_ch1",
  "gsm.extract_protocol_ch1",
  "gpl.technology"
)

sample_fields <- text_fields[text_fields != "gpl.technology"]

like_any_sql <- function(fields, terms) {
  clauses <- rep(
    paste0("lower(coalesce(", fields, ", '')) LIKE ?"),
    times = length(terms)
  )
  paste0("(", paste(clauses, collapse = " OR "), ")")
}

like_params <- function(terms, fields) {
  rep(paste0("%", str_to_lower(terms), "%"), each = length(fields))
}

has_one_mature_name <- function(text, mature_name) {
  lower_text <- str_to_lower(text)
  lower_name <- str_to_lower(mature_name)
  locations <- str_locate_all(lower_text, fixed(lower_name))[[1]]
  starts <- locations[, "start"]
  ends <- locations[, "end"]
  if (length(starts) == 0) {
    return(FALSE)
  }
  vapply(
    seq_along(starts),
    function(i) {
      before <- if (starts[i] == 1) "" else substr(lower_text, starts[i] - 1, starts[i] - 1)
      after <- if (ends[i] == nchar(lower_text)) "" else substr(lower_text, ends[i] + 1, ends[i] + 1)
      !str_detect(before, "[A-Za-z0-9]") && !str_detect(after, "[A-Za-z0-9]")
    },
    logical(1)
  ) |> any()
}

has_mature_name <- function(text, mature_name) {
  vapply(text, has_one_mature_name, mature_name = mature_name, logical(1))
}

classify_text <- function(text) {
  lower_text <- str_to_lower(text)
  kd <- knockdown_terms[str_detect(lower_text, fixed(knockdown_terms, ignore_case = TRUE))]
  ko <- knockout_terms[str_detect(lower_text, fixed(knockout_terms, ignore_case = TRUE))]
  oe <- overexpression_terms[str_detect(lower_text, fixed(overexpression_terms, ignore_case = TRUE))]
  class <- if (length(ko) > 0 && length(oe) == 0) {
    "knockout"
  } else if (length(kd) > 0 && length(oe) == 0) {
    "knockdown_or_inhibition"
  } else if ((length(kd) > 0 || length(ko) > 0) && length(oe) > 0) {
    "mixed_or_ambiguous"
  } else if (length(oe) > 0) {
    "overexpression"
  } else {
    "unknown"
  }
  list(class = class, terms = paste(c(kd, ko, oe), collapse = ","))
}

open_connection <- function() {
  dbConnect(SQLite(), db_path, flags = SQLITE_RO)
}

con <- open_connection()
on.exit(dbDisconnect(con), add = TRUE)

required_tables <- c("gse", "gsm", "gse_gsm", "gpl", "gse_gpl")
missing_tables <- setdiff(required_tables, dbListTables(con))
if (length(missing_tables) > 0) {
  stop("Database is missing required tables: ", paste(missing_tables, collapse = ", "))
}

sql <- paste0(
""
)

query_one <- function(search_terms, query_mirna = NA_character_, connection = con) {
  sql <- paste0(
    "
    SELECT DISTINCT
      gse.gse,
      gsm.gsm,
      gse.title AS gse_title,
      gsm.title AS gsm_title,
      gse.type AS gse_type,
      gse.submission_date,
      gse.last_update_date,
      gpl.gpl,
      gpl.title AS gpl_title,
      gpl.technology,
      gpl.organism AS gpl_organism,
      gsm.organism_ch1,
      gsm.source_name_ch1,
      gsm.characteristics_ch1,
      gsm.treatment_protocol_ch1,
      gse.summary,
      gse.overall_design
    FROM gse
    JOIN gse_gsm ON gse.gse = gse_gsm.gse
    JOIN gsm ON gse_gsm.gsm = gsm.gsm
    LEFT JOIN gse_gpl ON gse.gse = gse_gpl.gse
    LEFT JOIN gpl ON gse_gpl.gpl = gpl.gpl
    WHERE ",
    like_any_sql(text_fields, rna_seq_terms),
    " AND ",
    like_any_sql(sample_fields, search_terms),
    " AND ",
    like_any_sql(sample_fields, selected_perturbation_terms),
    if (length(date_years) > 0) {
      paste0(" AND substr(gse.submission_date, 1, 4) IN (", paste(rep("?", length(date_years)), collapse = ","), ")")
    } else {
      ""
    },
    if (nzchar(species_filter)) {
      " AND (
        lower(coalesce(gpl.organism, '')) LIKE ?
        OR lower(coalesce(gpl.title, '')) LIKE ?
        OR lower(coalesce(gsm.organism_ch1, '')) LIKE ?
        OR lower(coalesce(gsm.characteristics_ch1, '')) LIKE ?
        OR lower(coalesce(gsm.source_name_ch1, '')) LIKE ?
      )"
    } else {
      ""
    },
    "
    ORDER BY gse.gse, gsm.gsm"
  )

  params <- c(
    like_params(rna_seq_terms, text_fields),
    like_params(search_terms, sample_fields),
    like_params(selected_perturbation_terms, sample_fields),
    date_years,
    if (nzchar(species_filter)) rep(paste0("%", str_to_lower(species_filter), "%"), 5) else character()
  )

  hits <- dbGetQuery(connection, sql, params = params)
  if (nrow(hits) == 0) {
    return(list(clean = hits, ambiguous = hits))
  }

  combined <- apply(hits, 1, function(row) paste(row, collapse = " "))
  if (!is.na(query_mirna)) {
    keep <- has_mature_name(combined, query_mirna)
    hits <- hits[keep, , drop = FALSE]
    combined <- combined[keep]
  }
  if (nrow(hits) == 0) {
    return(list(clean = hits, ambiguous = hits))
  }

  classes <- lapply(combined, classify_text)
  hits$query_mirna <- query_mirna
  hits$perturbation_class <- vapply(classes, `[[`, character(1), "class")
  hits$matched_perturbation_terms <- vapply(classes, `[[`, character(1), "terms")
  hits$matched_targets <- vapply(
    combined,
    function(text) {
      paste(search_terms[str_detect(str_to_lower(text), fixed(str_to_lower(search_terms)))], collapse = ",")
    },
    character(1)
  )
  ambiguous <- hits[hits$perturbation_class == "mixed_or_ambiguous", , drop = FALSE]
  clean <- hits[hits$perturbation_class != "mixed_or_ambiguous", , drop = FALSE]
  list(clean = clean, ambiguous = ambiguous)
}

message("Searching ", length(mature_mirnas), " mature miRNAs")
search_mature_mirna <- function(mirna) {
  worker_con <- open_connection()
  on.exit(dbDisconnect(worker_con), add = TRUE)
  message("Searching ", mirna)
  query_one(mirna, query_mirna = mirna, connection = worker_con)
}

if (opts$cpus > 1) {
  cluster <- makeCluster(opts$cpus)
  on.exit(stopCluster(cluster), add = TRUE)
  clusterEvalQ(cluster, {
    library(DBI)
    library(RSQLite)
    library(stringr)
  })
  clusterExport(
    cluster,
    varlist = c(
      "db_path",
      "opts",
      "date_years",
      "species_filter",
      "rna_seq_terms",
      "knockdown_terms",
      "knockout_terms",
      "overexpression_terms",
      "selected_perturbation_terms",
      "text_fields",
      "sample_fields",
      "like_any_sql",
      "like_params",
      "has_one_mature_name",
      "has_mature_name",
      "classify_text",
      "open_connection",
      "query_one"
    ),
    envir = environment()
  )
  query_results <- parLapply(cluster, mature_mirnas, search_mature_mirna)
} else {
  query_results <- lapply(mature_mirnas, search_mature_mirna)
}
hits <- do.call(rbind, lapply(query_results, `[[`, "clean"))
ambiguous_hits <- do.call(rbind, lapply(query_results, `[[`, "ambiguous"))

format_hits <- function(rows) {
  if (nrow(rows) == 0 && !"query_mirna" %in% names(rows)) {
    return(rows)
  }
  rows <- rows[, c(
    "query_mirna",
    "gse",
    "gsm",
    "perturbation_class",
    "matched_targets",
    "matched_perturbation_terms",
    setdiff(names(rows), c(
      "query_mirna",
      "gse",
      "gsm",
      "perturbation_class",
      "matched_targets",
      "matched_perturbation_terms"
    ))
  )]

  rows[!duplicated(paste(rows$query_mirna, rows$gse, sep = "\t")), ]
}

hits <- format_hits(hits)
ambiguous_hits <- format_hits(ambiguous_hits)

write_tsv(hits, opts$output)
message("Wrote ", nrow(hits), " rows to ", opts$output)

if (nzchar(opts$ambiguous_output)) {
  write_tsv(ambiguous_hits, opts$ambiguous_output)
  message("Wrote ", nrow(ambiguous_hits), " ambiguous rows to ", opts$ambiguous_output)
}
