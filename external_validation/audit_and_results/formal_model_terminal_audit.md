# Formal model terminal audit

- Registry rows: 192; unique unit identifiers: 192.
- Model invoked: 192/192.
- Completed: 160/192; every completed row has stable code `OK` or `OK_NO_COMPARABLE_CREDIBLE_SETS`, a posterior artifact flag, and a non-empty RDS SHA-256.
- Failed: 32/192; every failure has stable code `E_SUSIE_NONCONVERGENCE`, no posterior artifact, and no RDS SHA-256.
- CAD: 96 completed (11 `OK`; 85 `OK_NO_COMPARABLE_CREDIBLE_SETS`).
- HF: 64 completed (4 `OK`; 60 `OK_NO_COMPARABLE_CREDIBLE_SETS`) and 32 failed under the locked 1,000-iteration ceiling.
- Each of the four tissue identifiers contributed 48 units and exactly eight HF nonconvergence failures.
- Runner, implementation-lock, and canonical-protocol hashes each have one unique value across all 192 rows.
- Audit scope: terminal workflow behavior only. Posterior contents were not inspected or interpreted.
