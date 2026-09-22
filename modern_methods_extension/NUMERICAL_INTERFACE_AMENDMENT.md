# Numerical-interface amendment during revision engineering

Before reading the new challenge performance table, an interface diagnostic
identified that float32 storage can make projected LD matrices slightly
indefinite (for example R01 minimum eigenvalue approximately -4.08e-8).
This is below the fail-closed gate's double-precision PSD tolerance. Neural
message passing uses float32, whereas scientific input validation should use
the source's float64 representation, before tensor conversion.

Retain and report the first float32 gate audit. Also run a separately labelled
source-float64 gate audit using the same 100 selected variants and the declared
simulation-only eigenvalue-floor policy applied to the archived source LD.
This is an explicit interface experiment, not an automatic repair after a gate
failure. Retain the same z vectors, frozen models, scores, thresholds and faults;
do not retrain or select examples. Supplied-LD permutations are applied in the
same fixed way to both representations. Compare the graph-input and gate-input
matrices numerically and record hashes. No claim of canonical input equivalence
is allowed without reporting their numerical difference.

The initial engineering harness also supplied a quantitative GWAS label to a
gate explicitly restricted to case-control GWAS. That metadata error correctly
produced E_GWAS_TRAIT_TYPE_INVALID for every fixture. Its outputs are retained
as engineering diagnostics, not scored as content-error detections. The corrected
contract probe supplies the gate's declared synthetic case-control metadata and
case fraction 0.2 to every fixture. This adapter test does not assert that the
Gaussian summary simulations are actual binary-phenotype GWAS observations.
