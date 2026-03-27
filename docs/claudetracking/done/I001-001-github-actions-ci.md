# I001-001: Add GitHub Actions CI pipeline

**Epic**: INFRA-001
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

Set up a GitHub Actions workflow that runs on push and pull request. Jobs:
`cargo fmt -- --check`, `cargo clippy -- -D warnings`, `cargo test`,
`cargo audit` (if available).

## Success Criteria

- [x] `.github/workflows/ci.yml` exists and is valid.
- [x] Workflow runs fmt, clippy, test, and audit.
- [x] Workflow triggers on push to main and on pull requests.
- [x] /update has been run after changes.
