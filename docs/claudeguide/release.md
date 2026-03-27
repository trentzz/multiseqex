# Release Configuration

Release configuration for the multiseqex crate.

## Registry

Published to [crates.io](https://crates.io/).

## Versioning

Follows the `v0.Y.Z` scheme. Bump the minor version for new features and the
patch version for bug fixes.

## Pre-publish Checks

Run all checks before publishing. Do not publish if any fail.

```sh
cargo fmt -- --check
cargo clippy -- -D warnings
cargo test
```

## Publish

Authentication is pre-configured via `cargo login`.

```sh
cargo publish
```

## Post-publish

Create a GitHub release after the crate is published.

```sh
gh release create vX.Y.Z --title "vX.Y.Z" --generate-notes
```

Replace `X.Y.Z` with the version being released.

## Package Exclusions

Only the following should be included in the published crate:

- Rust source files (`src/`)
- `Cargo.toml`
- `Cargo.lock`
- `README.md`
- `LICENSE`

Exclude everything else from the package. In particular, exclude:

- `docs/`
- `tests/`
- `docs/claudetracking/`

Use the `exclude` field in `Cargo.toml` to enforce this.

```toml
[package]
exclude = ["docs/", "tests/"]
```
