library(SNPRelate)

# --- Workspace Setup ---
cat("A window will open: please select the input RAW GDS file (vcf_data.gds)...\n")
raw_gds <- file.choose() 

work_dir <- dirname(raw_gds)
setwd(work_dir)
cat("Working directory set to:", work_dir, "\n")

genofile <- snpgdsOpen(raw_gds)
all_samples <- read.gdsn(index.gdsn(genofile, "sample.id"))

# Interactive sample removal
cat("\nA window will pop up. Select samples to REMOVE.\n")

selected_items <- select.list(
  all_samples, 
  multiple = TRUE, 
  title = "Select to REMOVE (Click OK/Cancel to keep all)",
  graphics = TRUE
)

if (length(selected_items) == 0) {
  cat("\nNo samples selected for removal. Task finished.\n")
  snpgdsClose(genofile)
} else {
  cat("\nRemoving samples:", paste(selected_items, collapse = ", "), "\n")
  
  good_samples <- setdiff(all_samples, selected_items)
  
  cat("\n--- Running LD-pruning for remaining samples ---\n")
  
  set.seed(1000)
  snpset_final <- snpgdsLDpruning(
    gdsobj = genofile,
    sample.id = good_samples, 
    autosome.only = FALSE,
    maf = 0.05,            
    missing.rate = 0.10,   
    method = "composite",
    ld.threshold = 0.8,
    verbose = TRUE,
    start.pos = "first"
  )
  pruned_snp_final <- unlist(unname(snpset_final))
  
  cat("\n--- Creating Clean GDS File ---\n")
  
  clean_gds <- "clean_vcf_data.gds"
  if (file.exists(clean_gds)) file.remove(clean_gds)
  
  snpgdsCreateGenoSet(
    src.fn = raw_gds, 
    dest.fn = clean_gds, 
    sample.id = good_samples,
    snp.id = pruned_snp_final 
  )
  
  snpgdsClose(genofile)
  
  cat("\nSUCCESS! Clean GDS file created with updated LD pruning.\n")
}
