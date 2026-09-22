# Revision interpretation ledger

1. The historical single-signal condition `summary_order_permuted` permutes
   trait-2 z values, not their identifiers, and did not invoke the gate. Its
   historical `contract_detectable=TRUE` flag was an incorrect manual annotation,
   not an observed gate result. Current code sets that annotation to FALSE.
   Observed PIP/PP.H4 outcomes are unchanged. An actual identifier-carrying row
   permutation is a different fault; those are checked directly in Additional 3.
2. Original historical output tables and nested hashes are preserved, not
   backdated or silently rewritten. Use SHA256SUMS_CURRENT.tsv for this revision.
3. Historical protocol and runtime-manifest timestamps have an unresolved
   chronology inconsistency disclosed in the manuscript. No public prospective
   preregistration or independently established result-naive timing is claimed.
4. The author reports HERMES conditions confirmed; the documentary record is
   pending. Third-party raw summary statistics/genotypes are not redistributed.
5. Review-only licence wording is superseded for project-owned files by the
   public Apache-2.0 / CC BY 4.0 scope; third-party rights remain unchanged.
