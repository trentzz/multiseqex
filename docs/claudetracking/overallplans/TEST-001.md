# TEST-001: Test Coverage Gaps

## Goal

Close remaining test coverage gaps found during the second audit.

## Motivation

The review cycle added good unit and integration tests, but several edge cases
remain untested: empty tables, malformed CSV rows, gzip rejection at CLI level,
combining multiple input sources, duplicate regions, and filename collisions.

## Scope

- Integration tests for untested paths.
- Edge-case unit tests.

## Tasks

- T001-001: Test empty table (headers only, no data rows)
- T001-002: Test malformed CSV rows (missing fields, bad values)
- T001-003: Integration test for gzip file rejection
- T001-004: Test combining --regions with --table
- T001-005: Test filename collision in --output-dir
- T001-006: Test --list with comments and blank lines
- T001-007: Test FAI build round-trip accuracy
- T001-008: Test --threads flag
