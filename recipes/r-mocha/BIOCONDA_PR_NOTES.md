# Bioconda submission: `r-mocha`

Use this text when opening a PR to [bioconda/bioconda-recipes](https://github.com/bioconda/bioconda-recipes).

## Package summary

- **Name:** `r-mocha`
- **Version:** 2.0.0
- **Upstream:** https://github.com/aifimmunology/MOCHA
- **License:** GPL-3.0-or-later (`LICENSE.md`)

## Source provenance

`MOCHA_2.0.0.tar.gz` is not present on active CRAN `src/contrib` (latest archived source is 1.1.0). This recipe pins an immutable GitHub commit tarball:

- **URL:** `https://github.com/aifimmunology/MOCHA/archive/93d513fdaf4d7c0f056d03ad7dc3f563e03713ef.tar.gz`
- **sha256:** `df3f539f88b86858e8963f1b7a04b2c6960bade26cb16fd5366ce6e5899a1c26`

When a tagged release tarball for 2.0.0 is published, prefer switching `source.url` to that release archive.

## Local validation

### Build

```bash
conda mambabuild recipes/r-mocha \
  --override-channels -c conda-forge -c bioconda -c defaults
```

**Result:** success (`exit 0`)

**Artifact:**

```
/home/enki/miniforge3/conda-bld/noarch/r-mocha-2.0.0-r45_0.tar.bz2
```

Recipe tests (during build): `library('MOCHA')`, `packageVersion('MOCHA')` → `2.0.0`.

### Clean-environment install smoke test

```bash
conda create -n r-mocha-smoke -y \
  --override-channels \
  -c file:///home/enki/miniforge3/conda-bld \
  -c conda-forge -c bioconda -c defaults \
  r-mocha

conda run -n r-mocha-smoke R -q -e "library('MOCHA'); print(packageVersion('MOCHA'))"
```

**Result:** `[1] ‘2.0.0’`

## Recipe notes for reviewers

- `noarch: generic` R package; `skip: true  # [win]`
- Runtime deps mirror `DESCRIPTION` `Imports` (Bioconductor + CRAN)
- Heavy `Suggests` (e.g. ArchR, BSgenome annotations) are intentionally not hard requirements
- Maintainer: `MPebworthEpana`

## Suggested PR checklist (upstream)

- [ ] Copy `recipes/r-mocha/meta.yaml` into bioconda-recipes fork
- [ ] Confirm maintainer GitHub handle
- [ ] Let Bioconda CI run (linux-64 / osx-64 as applicable)
- [ ] If solver fails, tighten `r-base` to one Bioconductor line in both `host` and `run`

## Post-merge install (users)

```bash
conda install -c conda-forge -c bioconda r-mocha
```
