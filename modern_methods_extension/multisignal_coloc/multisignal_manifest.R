list(schema_version = "1.0", seed_base = 20260920L, replicates = 50L, 
    variants = 100L, sample_sizes = c(2000L, 5000L, 20000L), 
    architectures = c("one_shared", "two_shared", "shared_plus_specific", 
    "distinct_only"), ld_policies = c("EUR_matched", "EAS_substituted", 
    "AFR_substituted"), data_generating_ld = "1000 Genomes EUR common-variant LD; nearest-PSD projection for simulation only", 
    consumers = c("susieR 0.14.2", "coloc 5.2.3"), model_settings = "SuSiE-RSS L=5, max_iter=1000, tol=1e-3; coloc priors p1=1e-4, p2=1e-4, p12=1e-5", 
    interpretation_boundary = "Synthetic robustness benchmark only. It compares consumer behavior across declared architectures, sample sizes and LD substitutions; it does not estimate empirical cardiovascular association, biological colocalisation, clinical validity or superiority over another software package.")
