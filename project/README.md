# Generated schema artifacts

The LinkML schema in `src/gocam/schema/gocam.yaml` is the source of truth for the GO-CAM data
model. Everything else in this directory is generated, disposable local output and is ignored by
Git.

Run `make gen-project` from the repository root to generate non-documentation schema artifacts here.
Run `make gendoc` to generate schema documentation separately. Versioned copies of this directory
are available as downloadable archives attached to GitHub Releases.
