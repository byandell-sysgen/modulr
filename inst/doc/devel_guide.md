# Developer Guide Blueprint, Implementation Plan & Walkthrough for `modulr`

This document records the developer guide prompt blueprint, implementation plan, and completion walkthrough for **`modulr`** (v0.3.1).

---

## 1. Actionable Prompt & Blueprint

### Purpose

Create a comprehensive developer guide infrastructure across R package repositories (`modulr`, `metabr`, `foundrHarmony`) detailing codebase architecture, WGCNA S3 object specifications, trait module data flow, and automated deployment via `pkgdown` and GitHub Actions.

### Blueprint Specifications

1. **Root Guide (`DEVELOPER.md`)**:
   - High-level package introduction, WGCNA ecosystem role, directory structure, function reference, developer quick start (`devtools::load_all()`), and release guidelines.

2. **Vignettes & Articles (`vignettes/devel_guide/`)**:
   - `index.Rmd`: Overview, ecosystem integration flowchart, and local development commands.
   - `modules.Rmd`: Technical reference detailing `wgcnaModules` objects, eigentraits (`out$eigen`), trait factors (`out$modules`), and `foundr` integration.
   - `data_flow.Rmd`: Data pipeline specifications covering matrix pivoting (`wgcna_pivot`), sample ID standardization (`wgcna_ID`), TOM dissimilarity math, and dynamic tree cut module merging.

3. **`pkgdown` Site Configuration (`_pkgdown.yml`)**:
   - Bootstrap 5 template.
   - Dynamic Mermaid JS v10 dynamic flowchart rendering header block.
   - Custom **Developer Guide** navbar dropdown menu linking to rendered HTML articles.
   - Categorized function reference sections.

4. **GitHub Actions Deployment (`.github/workflows/pkgdown.yaml`)**:
   - Automated workflow to build `pkgdown` site in CI/CD and deploy HTML assets directly to the **`gh-pages`** branch.
   - Configured `.gitignore` (`docs/`) and `.Rbuildignore`.

---

## 2. Implementation & Walkthrough

| File Path | Description |
| --- | --- |
| [`DEVELOPER.md`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/DEVELOPER.md) | Primary developer reference detailing WGCNA co-expression module algorithms, S3 object schemas, and `foundr` integration. |
| [`vignettes/devel_guide/index.Rmd`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/vignettes/devel_guide/index.Rmd) | High-level developer overview, package purpose, local quick start, and top-level end-to-end Mermaid flowchart. |
| [`vignettes/devel_guide/modules.Rmd`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/vignettes/devel_guide/modules.Rmd) | WGCNA trait modules, eigentraits (`out$eigen`), and kME factors (`out$modules`) specification. |
| [`vignettes/devel_guide/data_flow.Rmd`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/vignettes/devel_guide/data_flow.Rmd) | Data pipeline specifications covering matrix pivoting (`wgcna_pivot`), sample ID standardization (`wgcna_ID`), TOM dissimilarity math, and dynamic tree cut module merging. |
| [`_pkgdown.yml`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/_pkgdown.yml) | `pkgdown` site configuration with Bootstrap 5 and Mermaid JS rendering. |
| [`.github/workflows/pkgdown.yaml`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/.github/workflows/pkgdown.yaml) | Automated GitHub Actions workflow deploying site to `gh-pages` branch. |
| [`.gitignore`](file:///Users/brianyandell/Documents/Research/byandell-sysgen/modulr/.gitignore) | Added `docs/`. |
