# rerand Package Guidelines

## Project Scope

This project provides the R package for rerandomization designs and estimators. Maintain rigorous statistical reasoning, reproducible computational work, and consistency between R code, package documentation, tests, examples, and exported interfaces.

## Auxiliary Sources

The Codex project may include auxiliary sources configured when the project is created. These sources provide supporting context and may include related research materials, such as `Rerand-Review`.

Treat auxiliary sources as reference material unless the user explicitly asks for changes to them. When auxiliary sources conflict with files in this project, the current project is authoritative for the R package.

## Development Commands

Run these commands from the `rerand/` repository root:

```sh
Rscript -e 'devtools::document(".")'
Rscript -e 'Rcpp::compileAttributes(".")'
R CMD build .
R CMD check .
```

Use the generation commands before building and checking the package. After changing Roxygen comments or exported R code, run `devtools::document(".")`. After changing C++ code marked with `// [[Rcpp::export]]`, run `Rcpp::compileAttributes(".")`. Then run `R CMD build .` and `R CMD check .` for package-level validation.

## Generated Package Files

- Edit R source files under `R/`, C++ source files under `src/`, and documentation comments in those source files.
- `devtools::document(".")` regenerates `NAMESPACE` and files in `man/` from Roxygen comments.
- `Rcpp::compileAttributes(".")` regenerates `R/RcppExports.R` and `src/RcppExports.cpp` from exported C++ functions.
- Never manually edit generated `NAMESPACE`, `man/*.Rd`, `R/RcppExports.R`, or `src/RcppExports.cpp`. If a generated file appears incorrect, fix its source annotations or source code and regenerate it.

## Style and Testing

Use two-space indentation, snake_case for local variables, explicit package namespaces, and the established exported naming conventions. Keep Rcpp implementation details in `src/` aligned with the R-facing API. There is no formal `testthat` suite or coverage threshold.

- Use fixed seeds for smoke tests and numerical examples.
- Check both the R and Rcpp execution paths when a change affects computational behavior.
- Do not treat successful documentation generation as package validation; run `R CMD check .` for every substantive change.
