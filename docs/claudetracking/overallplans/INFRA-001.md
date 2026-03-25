# INFRA-001: CI and Packaging

## Goal

Set up continuous integration and improve packaging metadata.

## Motivation

No CI pipeline exists. Pull requests and pushes are not automatically tested.
The project also lacks a CHANGELOG and documented MSRV.

## Scope

- GitHub Actions CI (test, clippy, fmt, audit).
- CHANGELOG.md.
- Documented MSRV in Cargo.toml and README.

## Tasks

- I001-001: Add GitHub Actions CI pipeline
- I001-002: Add CHANGELOG.md
- I001-003: Document MSRV (Rust 1.85+)
