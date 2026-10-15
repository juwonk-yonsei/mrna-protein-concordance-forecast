<!--
B_revision_B15_claude.md — Claude's manuscript text for the final B stage, applied by the agent in B 지시서 15 §3.
  [OPTION Y/Z] choose by B 지시서 15 §2.4: Y if the CSV round trip is shown to explain the difference, Z if it does not.
Only one of H1/H2 is applied. If OPTION Z, nothing in this file changes the manuscript except H3.
-->

# H1. Supplementary Note 2 — the reproduction paragraph [OPTION Y only]

**Find** the sentence "We could not identify the cause." in the reproduction paragraph **and replace it with:**

The cause is the path the training rows took. The model frozen for RPPA and ProCan was fitted on the leave-one-cancer-out table as it stood in memory, while the APOLLO fit and the released model read the same table back from the tab-separated file it had just been written to. That round trip is not value-preserving at full precision [FILL: one clause naming what was lost, e.g. "because the file stores the feature columns to eight significant digits"], and gradient-boosted trees amplify such a difference: a feature value that crosses a histogram boundary changes a split, and the change propagates through the remaining trees. Refitting on the values as the file stores them reproduces the APOLLO predictions exactly, and refitting on [FILL: the recovered in-memory values, named] reproduces the frozen RPPA and ProCan predictions (maximum absolute differences [FILL] and [FILL]). Both fitted models are released.

<!-- src: results/stage15/roundtrip_check.tsv; results/stage15/roundtrip_refit.tsv -->

---

# H2. Methods, "Forecast models and baselines" — the two-fit paragraph [OPTION Y only]

**Find:** "The released model reproduces the APOLLO predictions exactly but not the RPPA and ProCan predictions (maximum absolute differences 0.031 and 0.044; Spearman correlations 0.9991 and 0.9992), and the cause was not identified (Supplementary Note 2)."

**Replace with:**

The two fits differ only in whether the leave-one-cancer-out table was taken from memory or read back from the file it had been written to; that difference in the stored precision of the features changes which splits the trees take (Supplementary Note 2). Each fit is reproduced exactly from its own input, and both models are released.

---

# H3. Supplementary Note 2 — one sentence to add in either case

**Append at the end of the reproduction paragraph** (after the sentence about the frozen predictions being released):

Code that fits a model on a table in memory and later reads the same table from disk should treat the two as different inputs.

This is the lesson a reader can use, and it holds whether or not the cause was confirmed. [OPTION Z: keep this sentence; it is the only change this file makes.]
