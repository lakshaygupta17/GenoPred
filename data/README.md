# Data directory

Large biomedical datasets are intentionally excluded from version control.

- `raw/`: immutable files downloaded from the named repository.
- `interim/`: harmonized coordinates, normalized tracks and cached tables.
- `processed/`: model-ready gene tables and graph artifacts.

Every downloaded file must be accompanied by a provenance record containing its source ID,
download URL, retrieval date, file size, SHA-256 checksum, access class and genome build.
Raw files are never edited in place.

