# Publishing

This repository is configured for GitHub Trusted Publishing with both TestPyPI
and PyPI via [`.github/workflows/publish.yml`](../.github/workflows/publish.yml).

## One-time setup on TestPyPI

If `gmxtopology` does not exist on TestPyPI yet, create a pending publisher:

- Go to <https://test.pypi.org/manage/account/publishing/>
- Add a GitHub publisher with:
  - Owner: `vojtechkostal`
  - Repository name: `gmxtopology`
  - Workflow filename: `publish.yml`
  - Environment name: `testpypi`

If the project already exists on TestPyPI, add the same publisher from the
project's **Manage** -> **Publishing** page instead.

## One-time setup on PyPI

If `gmxtopology` does not exist on PyPI yet, create a pending publisher:

- Go to <https://pypi.org/manage/account/publishing/>
- Add a GitHub publisher with:
  - Owner: `vojtechkostal`
  - Repository name: `gmxtopology`
  - Workflow filename: `publish.yml`
  - Environment name: `pypi`

If the project already exists on PyPI, add the same publisher from the
project's **Manage** -> **Publishing** page instead.

## GitHub environment setup

Create two GitHub environments in this repository:

- `testpypi`
- `pypi`

These environments are optional for Trusted Publishing itself, but recommended.
They let you require manual approvals or restrict who can trigger a production
publish.

## How to publish

### TestPyPI

Run the `Publish Python Package` workflow manually and choose
`repository=testpypi`.

### PyPI

Either:

- Run the same workflow manually and choose `repository=pypi`, or
- Create a GitHub release, which triggers the PyPI publish job automatically.

## Release flow

1. Prepare the version bump and changelog on a feature branch.
2. Run `python -m pytest`.
3. Run `python -m build` and `python -m twine check dist/*`.
4. Open a pull request and merge it after the `Checks` workflow passes.
5. Tag the release commit on `main`, for example:

   ```bash
   git switch main
   git pull --ff-only
   git tag -a v0.2.0 -m "gmxtopology 0.2.0"
   git push origin v0.2.0
   ```

6. Run the `Publish Python Package` workflow manually with
   `repository=testpypi`.
7. Verify installation from TestPyPI in a clean environment.
8. Create a GitHub release from the tag. Publishing that release triggers the
   production PyPI job automatically.

For the first release only, register the pending publishers and GitHub
environments described above before running the workflow.
