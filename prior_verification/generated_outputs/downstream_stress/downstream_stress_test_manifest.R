list(seed_rule = "20260831 + replicate index", replicates_per_condition = 100L, 
    variants = 100L, parallel_workers = 10L, nominal_sample_size = 5000, 
    simulated_truth = "one shared causal variant; EUR LD data-generating matrix", 
    consumers = c("susieR 0.14.2", "coloc.abf 5.2.3"), coloc_priors = "p1=1e-4, p2=1e-4, p12=1e-5; quantitative traits; sdY=1; synthetic MAF=0.25", 
    important_boundary = "The contract can stop identifier/order and declared ancestry/LD provenance mismatches. An internally corrupted effect sign with otherwise self-consistent metadata is intentionally treated as an undetectable content-level failure boundary, not as a prevented error.")
