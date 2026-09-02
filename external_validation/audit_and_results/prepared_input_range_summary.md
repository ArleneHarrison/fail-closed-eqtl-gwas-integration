# Prepared-input scalar range audit

This read-only descriptive audit streamed all 192 prepared `summary.csv` files
(389,322 unit rows) after formal preparation. Values are repeated across
factorial units where the same source records contribute to more than one
attempt; the table is therefore a range audit, not a count of independent
variants or participants.

- CAD retained the file-level denominator in every unit row: 34,541 cases and
  total N 296,525 (row case fraction 0.116485962397774).
- HF retained variant-level denominators. Across prepared unit rows, total N
  ranged from 60,286 to 1,685,710 (median 1,631,460), and case count ranged
  from 12,629 to 132,176 (median 122,599). The row-level case fraction ranged
  from 0.030754 to 0.258037 (median 0.076007).
- The downstream case-control model did not recompute its case fraction from
  these varying rows; it used the independently frozen study-level value
  0.0816774392949421. Its scalar N followed the separately declared within-unit
  median rule.
- In the 24 selected regional extracts, the preparation audit recorded zero
  invalid numeric/allele GWAS rows and zero exact GWAS duplicates; this does
  not assert that every row in the full deposited HERMES archive is valid.

No posterior or biological output was read to produce this audit.
