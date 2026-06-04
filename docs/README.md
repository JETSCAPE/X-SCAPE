# X-SCAPE documentation (local build)

This directory holds the source for the X-SCAPE documentation site:

- **MkDocs Material** — the narrative docs (framework, physics modules, guides).
  Config: [`../mkdocs.yml`](../mkdocs.yml); pages: the `*.md` files here.
- **Doxygen** — the C++ API reference. Config: [`Doxyfile`](Doxyfile).

> **CI deploy is currently disabled.** The GitHub Actions workflow is kept as
> `.github/workflows/docs.yml.disabled` (GitHub only runs `.yml`/`.yaml`
> files), so nothing is built or deployed automatically yet. The docs are built
> **locally** until it has been decided how to integrate this site with
> JETSCAPE. See "Re-enabling CI" below.

## Build locally

From the repository root:

```bash
./docs/build_docs.sh          # build into ./site
./docs/build_docs.sh serve    # build, then live-preview at http://localhost:8000
```

Or run the steps by hand:

```bash
python3 -m pip install -r docs/requirements.txt   # mkdocs-material, extensions
python3 -m mkdocs serve                            # live preview at :8000
# full static build (what CI would publish):
python3 -m mkdocs build --strict --site-dir site
doxygen docs/Doxyfile                              # adds the API under site/api/cpp
```

The built site lands in `./site/` (git-ignored — add `site/` to `.gitignore`
if you build in-tree). Open `site/index.html`, or use `mkdocs serve` for
auto-reload.

## Requirements

- Python 3 + `pip` (for `mkdocs-material`, see [`requirements.txt`](requirements.txt))
- `doxygen` + `graphviz` (the `dot` tool) for the API reference — optional; the
  build script skips the API step if `doxygen` is not installed.

## Re-enabling CI / GitHub Pages later

1. Rename `.github/workflows/docs.yml.disabled` back to
   `.github/workflows/docs.yml`.
2. In the repository: **Settings → Pages → Source → GitHub Actions**.
3. Push to `main` (or the configured branch). The site publishes at
   `https://<owner>.github.io/X-SCAPE/`.

The `site_url` and the absolute `/X-SCAPE/api/cpp/` links in the pages assume
the repository is served under the `X-SCAPE` path; adjust `mkdocs.yml` if the
docs are folded into a different JETSCAPE Pages site.
