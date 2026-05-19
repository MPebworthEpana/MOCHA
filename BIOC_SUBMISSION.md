# Opening the Bioconductor Contributions pull request

Use this checklist when `bioc/submission` passes `R CMD check` and `BiocCheck` on devel.

## Prerequisites

- [ ] Maintainer subscribed to [bioc-devel mailing list](https://stat.ethz.ch/mailman/listinfo/bioc-devel)
- [ ] `bioc/submission` pushed to `https://github.com/aifimmunology/MOCHA`
- [ ] `bioc-comments.md` updated with latest check output

## Contributions repository

1. Fork or use [Bioconductor/Contributions](https://github.com/Bioconductor/Contributions).
2. Open a **new issue** with the **New Package** template.
3. Provide:

| Field | Value |
|-------|--------|
| Package name | MOCHA |
| GitHub URL | https://github.com/aifimmunology/MOCHA |
| Branch | bioc/submission |
| Maintainer | Imran McGrath \<imran.mcgrath@alleninstitute.org\> |

4. Confirm the package builds on the [Single Package Builder](https://bioconductor.org/checkResults/) after the issue is opened.

## After acceptance

- Add remote: `git remote add upstream git@git.bioconductor.org:packages/MOCHA.git`
- Push devel: `git push upstream master:master` (from merged `main`)
- Create `RELEASE_3_xx` branch for each Bioconductor release cycle
