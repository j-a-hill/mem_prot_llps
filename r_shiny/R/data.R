LOCATION_CANONICAL <- list(
  list(kws = c("plasma membrane", "cell membrane"), name = "Plasma membrane"),
  list(kws = c("nucleus", "nucleoplasm", "nucleolus", "nuclear pore", "nuclear envelope"), name = "Nucleus"),
  list(kws = c("cytoplasm", "cytosol"), name = "Cytoplasm"),
  list(kws = c("endoplasmic reticulum"), name = "ER"),
  list(kws = c("golgi apparatus", "golgi"), name = "Golgi apparatus"),
  list(kws = c("mitochondri"), name = "Mitochondrion"),
  list(kws = c("lysosom"), name = "Lysosome"),
  list(kws = c("early endosome", "late endosome", "endosome"), name = "Endosome"),
  list(kws = c("peroxisom"), name = "Peroxisome"),
  list(kws = c("synapse", "presynaptic", "postsynaptic", "synaptic vesicle"), name = "Synapse"),
  list(kws = c("cell junction", "tight junction", "adherens junction", "gap junction"), name = "Cell junction"),
  list(kws = c("secreted", "extracellular space", "extracellular region"), name = "Extracellular"),
  list(kws = c("cell surface"), name = "Cell surface"),
  list(kws = c("stress granule"), name = "Stress granule"),
  list(kws = c("p-body", "processing body"), name = "P-body"),
  list(kws = c("lipid droplet"), name = "Lipid droplet"),
  list(kws = c("autophagosom", "phagosom"), name = "Autophagosome"),
  list(kws = c("vesicle"), name = "Vesicle"),
  list(kws = c("centrosom", "centriole", "spindle pole"), name = "Centrosome"),
  list(kws = c("axon", "dendrite", "neurite", "growth cone"), name = "Neuronal process"),
  list(kws = c("cilium", "flagellum", "basal body"), name = "Cilium"),
  list(kws = c("chromatin", "chromosome", "kinetochore"), name = "Chromatin/Chromosome"),
  list(kws = c("focal adhesion", "focal contact"), name = "Focal adhesion"),
  list(kws = c("sarcomere", "z disc", "myofibril"), name = "Sarcomere"),
  list(kws = c("exosom"), name = "Exosome")
)

FUNCTION_CATEGORIES <- list(
  "Ion channel" = c(
    "\\bion\\s*channel\\b", "\\bsodium\\s+channel\\b", "\\bpotassium\\s+channel\\b",
    "\\bcalcium\\s+channel\\b", "\\bchloride\\s+channel\\b", "\\bchannel\\s+protein\\b",
    "\\bvoltage.gated\\b", "\\bligand.gated\\b", "\\bcation\\s+channel\\b",
    "\\banion\\s+channel\\b", "\\baquaporin\\b", "\\bgap\\s+junction\\b",
    "\\bconnexin\\b", "\\bpannexin\\b"
  ),
  "GPCR" = c(
    "\\bg\\s+protein.coupled\\s+receptor\\b", "\\bgpcr\\b",
    "\\b7.transmembrane\\b", "\\b7tm\\b", "\\bmetabotropic\\s+receptor\\b",
    "\\bserpentine\\s+receptor\\b"
  ),
  "Receptor tyrosine kinase" = c(
    "\\breceptor\\s+tyrosine\\s+kinase\\b", "\\brtk\\b",
    "\\btyrosine.protein\\s+kinase.*receptor\\b",
    "\\begf\\s+receptor\\b", "\\bfgf\\s+receptor\\b", "\\bpdgf\\s+receptor\\b",
    "\\bvegf\\s+receptor\\b", "\\binsulin\\s+receptor\\b"
  ),
  "Transporter" = c(
    "\\btransporter\\b", "\\babc\\s+transporter\\b", "\\bsolute\\s+carrier\\b",
    "\\bslc\\b", "\\bsymporter\\b", "\\bantiporter\\b", "\\buniporter\\b",
    "\\bpermease\\b", "\\bcarrier\\s+protein\\b", "\\bimporter\\b", "\\bexporter\\b"
  ),
  "Kinase" = c(
    "\\bkinase\\b", "protein\\s+kinase", "\\bphosphorylates\\b",
    "\\bcdk\\b", "\\bcyclin.dependent\\s+kinase\\b", "\\bmapk\\b"
  ),
  "Phosphatase" = c(
    "\\bphosphatase\\b", "\\bprotein\\s+phosphatase\\b", "\\bdephosphorylates\\b", "\\bpten\\b"
  ),
  "Protease" = c(
    "\\bprotease\\b", "\\bpeptidase\\b", "\\bcaspase\\b",
    "\\bproteolytic\\b", "\\bmetalloprotease\\b", "\\bserine\\s+protease\\b",
    "\\bcysteine\\s+protease\\b", "\\bprotein\\s+cleavage\\b"
  ),
  "Ligase" = c(
    "\\bligase\\b", "\\bubiquitin\\s+ligase\\b", "\\be3\\s+ligase\\b",
    "\\bubiquitination\\b", "\\bsumo\\s+ligase\\b"
  ),
  "Synthetase" = c(
    "\\bsynthetase\\b", "\\bsynthase\\b", "\\baminoacyl.trna\\s+synthetase\\b"
  ),
  "Transferase" = c(
    "\\btransferase\\b", "\\bmethyltransferase\\b", "\\bacetyltransferase\\b",
    "\\bglycosyltransferase\\b", "\\bphosphotransferase\\b"
  ),
  "Oxidoreductase" = c(
    "\\boxidoreductase\\b", "\\bdehydrogenase\\b", "\\boxidase\\b",
    "\\breductase\\b", "\\bperoxidase\\b", "\\bcytochrome\\b"
  ),
  "Hydrolase" = c(
    "\\bhydrolase\\b", "\\besterase\\b", "\\blipase\\b", "\\batpase\\b", "\\bhydrolytic\\b"
  ),
  "Receptor" = c("\\breceptor\\b", "\\breceptor\\s+activity\\b"),
  "Nuclear receptor" = c(
    "\\bnuclear\\s+receptor\\b", "\\bsteroid.*receptor\\b",
    "\\bestrogen\\s+receptor\\b", "\\bandrogen\\s+receptor\\b",
    "\\bretinoic\\s+acid\\s+receptor\\b", "\\bvitamin\\s+d\\s+receptor\\b"
  ),
  "Transcription factor" = c(
    "\\btranscription\\s+factor\\b", "\\btranscriptional.*regulator\\b",
    "\\btranscriptional.*activator\\b", "\\btranscriptional.*repressor\\b",
    "\\bdna.binding.*transcription\\b"
  ),
  "GTPase" = c(
    "\\bgtpase\\b", "\\bras\\s+protein\\b", "\\brho\\s+protein\\b",
    "\\brab\\s+protein\\b", "\\bguanine\\s+nucleotide.binding\\b"
  ),
  "Structural" = c(
    "\\bcytoskeleton\\b", "\\bactin\\b", "\\btubulin\\b",
    "\\bintermediate\\s+filament\\b", "\\bcollagen\\b", "\\blaminin\\b",
    "\\bfibronectin\\b", "\\bspectrin\\b", "\\bscaffold\\b"
  ),
  "Adhesion" = c(
    "\\badhesion\\b", "\\bcadherin\\b", "\\bintegrin\\b",
    "\\bselectin\\b", "\\bcell.cell.*junction\\b", "\\bcell.*adhesion\\b"
  ),
  "Chaperone" = c(
    "\\bchaperone\\b", "\\bheat\\s+shock\\s+protein\\b", "\\bhsp\\b",
    "\\bprotein\\s+folding\\b", "\\bchaperonin\\b", "\\bco.chaperone\\b"
  ),
  "RNA processing" = c(
    "\\brna.binding\\b", "\\brna\\s+binding\\b", "\\brna\\s+processing\\b",
    "\\brna\\s+splicing\\b", "\\bpre.mrna.*splicing\\b", "\\bsplicing\\s+factor\\b",
    "\\bribosom\\w+\\b", "\\brna\\s+helicase\\b", "\\btranslation.*factor\\b",
    "\\btranslational.*regulator\\b"
  ),
  "DNA repair" = c(
    "\\bdna\\s+repair\\b", "\\bdna\\s+damage\\b", "\\bdna\\s+replication\\b",
    "\\bhomologous\\s+recombination\\b", "\\bgenome.*stability\\b",
    "\\bdna\\s+damage\\s+response\\b", "\\bdouble.strand\\s+break\\b"
  ),
  "Chromatin remodeling" = c(
    "\\bchromatin\\s+remodel\\w*\\b", "\\bhistone.*methyltransferase\\b",
    "\\bhistone.*acetyltransferase\\b", "\\bhistone.*deacetylase\\b",
    "\\bhistone.*demethylase\\b", "\\bepigenetic\\b", "\\bdna\\s+methyltransferase\\b",
    "\\bnucleosome.*remodeling\\b"
  )
)


parse_location <- function(raw) {
  if (is.na(raw) || !nzchar(trimws(raw))) return(character(0))
  text <- gsub("SUBCELLULAR LOCATION:", "", raw, ignore.case = TRUE)
  text <- gsub("\\{[^}]*\\}", "", text)
  text <- gsub("\\[[^\\]]*\\]", "", text)
  text <- gsub("\\([^)]*\\)", "", text)
  text <- gsub("Note=.*", "", text, ignore.case = TRUE)
  text <- tolower(text)
  found <- character(0)
  for (entry in LOCATION_CANONICAL) {
    for (kw in entry$kws) {
      if (grepl(kw, text, fixed = TRUE)) {
        if (!(entry$name %in% found)) found <- c(found, entry$name)
        break
      }
    }
  }
  found
}


parse_functions <- function(func_str, name_str = NA_character_) {
  parts <- character(0)
  if (!is.na(func_str) && nzchar(func_str)) parts <- c(parts, func_str)
  if (!is.na(name_str) && nzchar(name_str)) parts <- c(parts, name_str)
  if (length(parts) == 0L) return(character(0))
  combined <- gsub("\\{[^}]*\\}", "", paste(parts, collapse = " "))
  combined <- tolower(combined)
  found <- character(0)
  for (cat_name in names(FUNCTION_CATEGORIES)) {
    for (pat in FUNCTION_CATEGORIES[[cat_name]]) {
      if (grepl(pat, combined, perl = TRUE, ignore.case = TRUE)) {
        found <- c(found, cat_name)
        break
      }
    }
  }
  found
}


count_tm_domains <- function(domain_str) {
  if (is.na(domain_str) || !nzchar(domain_str)) return(0L)
  m <- gregexpr("\\d+\\.\\.\\d+", domain_str, perl = TRUE)[[1L]]
  if (m[[1L]] == -1L) 0L else length(m)
}


enrich_df <- function(df) {
  if ("Subcellular location [CC]" %in% names(df)) {
    df$location_categories <- lapply(df[["Subcellular location [CC]"]], parse_location)
  } else {
    df$location_categories <- vector("list", nrow(df))
  }

  func_col <- if ("Function [CC]"  %in% names(df)) df[["Function [CC]"]]  else rep(NA_character_, nrow(df))
  name_col <- if ("Protein names"  %in% names(df)) df[["Protein names"]]  else rep(NA_character_, nrow(df))
  df$function_categories <- mapply(parse_functions, func_col, name_col, SIMPLIFY = FALSE)

  if ("p(LLPS)" %in% names(df)) {
    pllps <- df[["p(LLPS)"]]
    df$pLLPS_class <- factor(
      ifelse(pllps >= 0.7, "High", ifelse(pllps >= 0.4, "Medium", "Low")),
      levels = c("High", "Medium", "Low")
    )
  }

  if (!"TMD_count" %in% names(df)) {
    tmd <- if ("Transmembrane" %in% names(df)) sapply(df$Transmembrane, count_tm_domains) else integer(nrow(df))
    imd <- if ("Intramembrane" %in% names(df)) sapply(df$Intramembrane, count_tm_domains) else integer(nrow(df))
    df$TMD_count <- tmd + imd
  }

  df
}


load_data <- function(path) {
  ext <- tolower(tools::file_ext(path))
  df <- if (ext == "xlsx") {
    as.data.frame(readxl::read_excel(path))
  } else {
    as.data.frame(readr::read_csv(path, show_col_types = FALSE))
  }
  enrich_df(df)
}


load_default <- function() {
  candidates <- c(
    "data/full_dataset.csv",
    "data/sample_data.csv",
    "../results/full_dataset.csv",
    "../dashboard/full_dataset.csv"
  )
  path <- Find(file.exists, candidates)
  if (is.null(path)) stop("Could not find full_dataset.csv. Place it in r_shiny/data/")
  load_data(path)
}
