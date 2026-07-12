# Zenodo release guide — MontuPython

This document describes how to archive a MontuPython release on [Zenodo](https://zenodo.org/) and obtain a DOI. The repository already includes:

| File | Purpose |
|------|---------|
| `.zenodo.json` | Metadata pre-filled for each GitHub release (Zenodo reads this **first**) |
| `CITATION.cff` | Citation metadata for GitHub and other tools (update DOI after publishing) |

> **Note:** If both files exist, Zenodo uses **only** `.zenodo.json` when archiving from GitHub. Keep both in sync manually.

---

## One-time setup (first release only)

### 1. Create a Zenodo account

1. Go to [https://zenodo.org/](https://zenodo.org/) and sign in (preferably with your GitHub account).

### 2. Connect GitHub to Zenodo

1. Zenodo → **Account** (top right) → **GitHub** → **Connect**.
2. Authorize Zenodo on GitHub.
3. Under **Enable repositories**, switch on **`seap-udea/MontuPython`**.

Zenodo will now listen for new **GitHub Releases** on that repository.

### 3. Commit metadata files

Ensure `.zenodo.json` and `CITATION.cff` are on the default branch (`main`):

```bash
git add .zenodo.json CITATION.cff docs/zenodo-release.md
git commit -m "Add Zenodo and citation metadata for v0.10.0"
git push origin main
```

---

## Pre-release checklist

Before creating the GitHub release:

- [ ] Version bumped in `setup.py`, `montu/version.py`, and `.zenodo.json` (`version` + `title`)
- [ ] `CITATION.cff` updated (`version`, `date-released`)
- [ ] `WHATSNEW.md` updated for the new version
- [ ] `README.ipynb` exported: `make readme` (updates `README.md` for PyPI)
- [ ] Tests pass: `make test`
- [ ] Working tree clean: `git status`
- [ ] (Optional) PyPI upload: `make release RELMODE=release VERSION=x.y.z`

Current target version: **0.10.0**

---

## Create the GitHub release (triggers Zenodo)

Zenodo archives **GitHub Releases**, not plain git tags. Use the web UI or `gh`:

### Option A — GitHub web UI

1. Open [https://github.com/seap-udea/MontuPython/releases/new](https://github.com/seap-udea/MontuPython/releases/new)
2. **Choose a tag:** `v0.10.0` (create new tag from `main`)
3. **Release title:** `v0.10.0`
4. **Description:** paste the `## Version 0.10.*` section from `WHATSNEW.md`
5. Click **Publish release**

### Option B — GitHub CLI

```bash
git tag -a v0.10.0 -m "MontuPython v0.10.0"
git push origin v0.10.0

gh release create v0.10.0 \
  --title "v0.10.0" \
  --notes-file WHATSNEW.md
```

Within a few minutes Zenodo creates a new **draft** deposit linked to this release.

---

## Publish on Zenodo

1. Go to [https://zenodo.org/](https://zenodo.org/) → **Upload** → **Drafts** (or check your email).
2. Open the new draft for `seap-udea/MontuPython` / `v0.10.0`.
3. Verify metadata (pre-filled from `.zenodo.json`):
   - Title, creators, keywords, license (MIT)
   - Related links (GitHub, PyPI, docs, web app)
4. Click **Publish**.
5. Copy the DOI, e.g. `10.5281/zenodo.1234567`.

### Reserve DOI before publishing (optional)

On the draft page you can **Reserve a DOI** and add it to `CITATION.cff` and `README.ipynb` **before** clicking Publish.

---

## After publishing — update the repository

### 1. `CITATION.cff`

Uncomment and fill the `doi` and `preferred-citation` block:

```yaml
doi: 10.5281/zenodo.XXXXXXX
```

### 2. `README.ipynb` badge

Replace the placeholder DOI badge (cell 1):

```markdown
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
```

Then run `make readme` to refresh `README.md`.

### 3. `.zenodo.json` for the *next* release

Update `version`, `title`, and `description` before the next GitHub release.

### 4. Commit post-release updates

```bash
git add CITATION.cff README.ipynb README.md
git commit -m "Add Zenodo DOI for v0.10.0"
git push origin main
```

---

## Subsequent releases

Each new version needs:

1. Bump versions (see checklist above).
2. Update `.zenodo.json` and `CITATION.cff`.
3. Create a new GitHub release with tag `vX.Y.Z`.
4. Review and **Publish** the new Zenodo draft.
5. Zenodo issues a **versioned DOI** under the same concept DOI (e.g. `10.5281/zenodo.1234567` → `10.5281/zenodo.1234568`).

---

## Troubleshooting

| Problem | Likely cause | Fix |
|---------|--------------|-----|
| No draft on Zenodo | GitHub integration disabled or release not published | Re-enable repo in Zenodo → GitHub settings |
| Metadata wrong | Stale `.zenodo.json` on the tagged commit | Ensure metadata is committed **before** tagging |
| Zenodo import failed | Invalid `.zenodo.json` | Validate JSON; check `upload_type`, `license`, creator names |
| Badge still says “pending” | DOI not yet in README | Update badge URL after publish |

---

## Suggested citation (after DOI is assigned)

> Zuluaga, J. I. (2026). *MontuPython v0.10.0: astronomical ephemerides for the ancient world* (Version 0.10.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX
