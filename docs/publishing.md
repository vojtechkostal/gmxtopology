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

## Recommended first release flow

1. Register the pending publisher on TestPyPI.
2. Run the workflow manually with `repository=testpypi`.
3. Verify installation from TestPyPI.
4. Register the pending publisher on PyPI.
5. Publish to PyPI manually or by creating a GitHub release.
