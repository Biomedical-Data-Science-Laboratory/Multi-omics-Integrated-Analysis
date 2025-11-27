library(MOFA2)

# rna_mirna_protein_RNAseq_mofa: 481 samples * 16016 features
# filter_rna_mirna_protein_mirna_mofa: 481 samples * 195 features
# rna_mirna_protein_protein_mofa: 481 samples * 463 features

rna_mirna_protein_mofa_list <- list(
  RNAseq = rna_mirna_protein_RNAseq_mofa,
  miRNA = filter_rna_mirna_protein_mirna_mofa,
  Protein = rna_mirna_protein_protein_mofa)

# Create the MOFA object
rna_mirna_protein_mofa_obj <- create_mofa(rna_mirna_protein_mofa_list)

# Define options
rna_mirna_protein_data_opts  <- get_default_data_options(rna_mirna_protein_mofa_obj)
rna_mirna_protein_model_opts <- get_default_model_options(rna_mirna_protein_mofa_obj)
rna_mirna_protein_train_opts <- get_default_training_options(rna_mirna_protein_mofa_obj)

rna_mirna_protein_model_opts$num_factors <- 15
## default: 15

# Build and train the MOFA object
rna_mirna_protein_mofa_obj <- prepare_mofa(
  object         = rna_mirna_protein_mofa_obj,
  data_options     = rna_mirna_protein_data_opts,
  model_options    = rna_mirna_protein_model_opts,
  training_options = rna_mirna_protein_train_opts
)

rna_mirna_protein_outfile <- file.path(getwd(),"rna_mirna_protein_model.hdf5")
rna_mirna_protein_mofa_trained <- run_mofa(rna_mirna_protein_mofa_obj, rna_mirna_protein_outfile, use_basilisk = TRUE)

# Visualization

# Data overview
plot_data_overview(rna_mirna_protein_mofa_trained)
# Distribution heatmap
plot_variance_explained(rna_mirna_protein_mofa_trained)
# Bar plot
plot_variance_explained(rna_mirna_protein_mofa_trained, x="group", y="factor", plot_total = T)
