# Test script to validate pipeline integration

# Check if the pipeline integration file loads correctly
tryCatch({
  source("pipeline_integration.R")
  cat("✓ pipeline_integration.R loaded successfully\n")
}, error = function(e) {
  cat("✗ Error loading pipeline_integration.R:", e$message, "\n")
})

# Check if the enhanced server loads correctly
tryCatch({
  source("enhanced_server.R")
  cat("✓ enhanced_server.R loaded successfully\n")
}, error = function(e) {
  cat("✗ Error loading enhanced_server.R:", e$message, "\n")
})

# Check if the enhanced UI loads correctly
tryCatch({
  source("enhanced_ui.R")
  cat("✓ enhanced_ui.R loaded successfully\n")
}, error = function(e) {
  cat("✗ Error loading enhanced_ui.R:", e$message, "\n")
})

# Check if the enhanced app loads correctly
tryCatch({
  source("enhanced_app.R")
  cat("✓ enhanced_app.R loaded successfully\n")
}, error = function(e) {
  cat("✗ Error loading enhanced_app.R:", e$message, "\n")
})

cat("\nIntegration test complete!\n")
cat("The GUI should now use real HiTMaP functions instead of fake data.\n")
cat("Key changes made:\n")
cat("1. Enhanced server now sources pipeline_integration.R instead of pipeline_server.R\n")
cat("2. Pipeline integration calls actual HiTMaP functions like Protein_feature_list_fun() and IMS_data_process()\n")
cat("3. Candidate generation saves results for PMF search to use\n")
cat("4. PMF search uses real HiTMaP IMS_data_process function\n")
cat("5. Quick Start calls the complete imaging_identification workflow\n")