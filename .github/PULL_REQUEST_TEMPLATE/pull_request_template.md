Many thanks to contributing to nf-core/eager!

Please fill in the appropriate checklist below (delete whatever is not relevant). These are the most common things requested on pull requests (PRs).

## PR checklist

 - [ ] This comment contains a description of changes (with reason).
 - [ ] If you've fixed a bug or added code that should be tested, add tests!
   - [ ] If you've added a new tool - add to the software_versions process and a regex to `scrape_software_versions.py`
   - [ ] If necessary, also make a PR on the [nf-core/eager branch on the nf-core/test-datasets repo]( https://github.com/nf-core/test-datasets/pull/new/nf-core/eager).
 - [ ] Make sure your code lints (`nf-core lint .`).
 - [ ] Ensure the test suite passes (`nextflow run . -profile test,docker`).
 - [ ] Usage Documentation in `docs/usage.md` is updated.
 - [ ] Output Documentation in `docs/output.md` is updated.
 - [ ] `CHANGELOG.md` is updated.
 - [ ] `README.md` is updated (including new tool citations and authors/contributors).

**Learn more about contributing:** https://github.com/nf-core/eager/tree/master/.github/CONTRIBUTING.md
