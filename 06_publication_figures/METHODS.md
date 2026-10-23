# Figure calculations

All panels load tab-separated proximity results and retain only rows marked `analyzed`. The plotting scripts use the Agg backend, Times New Roman, vector PDF output, and PNG output at 600 dpi.

| Script | Calculation |
| --- | --- |
| `figure_4_12_rank_stability.py` | Inner-join the two rankings by exact drug name; correlate the pandas ranks of their rank values and plot paired ranks. |
| `figure_4_12_topn_stability.py` | For N = 50, 100, 154, 230, 307, and 500, form drug sets with rank <= N; plot intersection size divided by N and the intersection divided by union size. |
| `network_proximity_ranking_top20.py` | Sort by ascending Z-score and then drug name; assign consecutive plot ranks and highlight the first 20 drugs. |
| `figure_4_13_roc_clinical_validation.py` | Classify drugs using the supplied clinical-reference set; use negative Z-score as the prediction score; aggregate tied scores at their endpoints and integrate the ROC curve by the trapezoid rule. |
| `figure_4_13_roc_clinical_validation_grid.py` | Use the same ROC calculation with an added light grid. |
| `figure_4_13_clinical_evidence_z_distribution.py` | Plot Z-score boxplots for reference-positive and remaining drugs without boxplot outliers; overlay all points using horizontal Gaussian jitter from NumPy's default generator seeded with 42. |
| `figure_4_14_target_count_vs_z.py` | Plot `n_targets_used` against Z-score and correlate the pandas ranks of both variables. |
| `figure_4_14_exact_profile_redundancy.py` | Read workbook columns A:C; collect unique non-missing STRING IDs per normalized drug name; group analyzed drugs by their sorted target-ID tuples; count profiles shared by at least two drugs and plot their frequencies on a logarithmic y-axis. |
| `figure_4_15_screening_depth.py` | Sort by supplied rank; select the top 10%, 15%, and 20%, with counts rounded using floor(n*p/100 + 0.5); divide recovered reference-positive rows by all reference-positive analyzed rows. |

Clinical labels use an exact match after trimming and lowercasing the drug name against `CT_POSITIVE` in `figure_scripts/common.py`. All unmatched drugs form the background class; absence from this list does not establish negative clinical evidence. The 24 supplied reference names are unchanged. No new literature review or clinical curation was performed.

Exact-profile redundancy uses the full supplied mapping workbook, rather than the threshold-filtered targets used in proximity calculations. Drug names are trimmed and lowercased. Missing workbook drug names and missing STRING IDs are skipped. A drug without mapped targets is assigned an empty profile by the supplied implementation. These rules are preserved.

Input limitations remain relevant when adapting the scripts: ranking files should have one analyzed row per drug, ROC inputs need both label classes and finite Z-scores, and Top-N comparisons require a non-empty union. No input-cleaning or statistical rule was changed during publication preparation.
