# Manifests

| File | Contents |
|---|---|
| `SHA256SUMS` | SHA-256 of every file in this repository (`sha256sum -c manifests/SHA256SUMS`) |
| `files.csv` | Repository path, original project path, bytes, SHA-256, and whether the file appears in the 766-file freeze manifest of the validation and still matches it |
| `stage24_revision_manifest.csv`, `stage24_freeze_diff.json` | Hash records of the final manuscript revision (Stage 24) |
| `integrity/` | Integrity-check records: the end of validation (all 766 frozen files unchanged), and the final pre-submission check, where all 34 changes are manuscript-side |
| `historical/` | Artifact manifests of earlier stages |

The freeze manifest itself is at `../numerical_validation/adversarial_overnight/freeze_manifest.json`, with its preregistration hash in `prereg_hash.json`. Hashes establish file identity, not mathematical validity.
