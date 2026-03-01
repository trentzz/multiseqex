# Feature: Change to Column Names

## Problem

Currently, `--table` and `--sv-table` determine parsing behavior based on the
**number of columns** (2/3/4 for `--table`, 4/5/6/7 for `--sv-table`). This is
fragile and inflexible:

- Column order is implicit and must be memorized
- Extra metadata columns (e.g. GENE) cause parsing failures
- Adding a NAME column changes the "mode" and shifts all other column positions
- Users cannot reorder columns in their files

## Solution

Switch to **header-name-based column lookup**. The tool reads the CSV/TSV header
row, finds required columns by name, and ignores any extra columns. Column order
no longer matters.

## Expected Column Names

### `--table` (regular regions)

| Column   | Required?                        | Description                      |
|----------|----------------------------------|----------------------------------|
| `CHROM`  | Yes                              | Chromosome / contig name         |
| `START`  | Yes (range mode)                 | 1-based inclusive start position |
| `END`    | Yes (range mode)                 | 1-based inclusive end position   |
| `POS`    | Yes (position mode, needs flank) | Single coordinate position       |
| `NAME`   | No                               | Region label for output naming   |

**Mode detection:**
- If `START` + `END` columns exist -> range mode (flank optional, adds padding)
- If `POS` column exists (no `START`/`END`) -> position mode (requires `--flank`)
- If both `POS` and `START`/`END` exist -> error (ambiguous)

### `--sv-table` (structural variant regions)

| Column        | Required?                        | Description                      |
|---------------|----------------------------------|----------------------------------|
| `CHROM_LEFT`  | Yes                              | Left breakpoint chromosome       |
| `START_LEFT`  | Yes (range mode)                 | Left breakpoint start            |
| `END_LEFT`    | Yes (range mode)                 | Left breakpoint end              |
| `POS_LEFT`    | Yes (position mode, needs flank) | Left breakpoint position         |
| `CHROM_RIGHT` | Yes                              | Right breakpoint chromosome      |
| `START_RIGHT` | Yes (range mode)                 | Right breakpoint start           |
| `END_RIGHT`   | Yes (range mode)                 | Right breakpoint end             |
| `POS_RIGHT`   | Yes (position mode, needs flank) | Right breakpoint position        |
| `NAME`        | No                               | SV identifier for output naming  |

**Mode detection:**
- If `START_LEFT` + `END_LEFT` + `START_RIGHT` + `END_RIGHT` exist -> range mode
- If `POS_LEFT` + `POS_RIGHT` exist -> position mode (requires `--flank`)

## Implementation TODO

- [x] Create featuredocs directory and plan document
- [x] Rewrite `parse_regions_table` to use header-name lookup
  - Read headers, build name->index map
  - Detect mode from which columns are present
  - Validate required columns exist, give clear errors if missing
  - Look up values by column name instead of hardcoded index
  - Ignore extra columns gracefully
- [x] Rewrite `parse_regions_sv_table` to use header-name lookup
  - Same approach as above for SV-specific column names
- [x] Update CLI argument help text
  - `--table` description: list expected column names
  - `--sv-table` description: list expected column names
  - `--flank` description: clarify when it's required
- [x] Update README.md
  - Document expected column names with tables
  - Add example CSV/TSV snippets
  - Update usage block
- [x] Build, test, commit, and create PR into `development`

## Design Decisions

- **Case-insensitive matching**: Header names are matched case-insensitively
  (e.g. `chrom`, `Chrom`, `CHROM` all work)
- **Extra columns ignored**: Any columns not in the expected set are silently
  skipped, so users can keep metadata like GENE, STRAND, etc.
- **Clear error messages**: If required columns are missing, the error lists
  exactly which columns were expected and which were found
- **No backward compatibility shims**: The old column-count behavior is fully
  replaced since this is a new feature branch
