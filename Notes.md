# How to publish

## Pypi

Create update:

1) Change version number in [pyproject.toml](pyproject.toml)
2) Create new package and upload:

```bash
SCOARY_VERSION="?.?.?"
# build: will create files in dist/
poetry build
# test: install .whl file
pip install -U dist/scoary_2-$SCOARY_VERSION-py3-none-any.whl
# upload
poetry publish
```

## Docker / Podman

If you use docker, simply replace each `podman` with `docker`.

```shell
podman build --build-arg SCOARY_VERSION=$SCOARY_VERSION --tag troder/scoary-2 .
```

Publish docker image:

```shell
# podman login docker.io --get-login
# podman login docker.io
podman tag troder/scoary-2 troder/scoary-2:$SCOARY_VERSION
podman push troder/scoary-2:$SCOARY_VERSION

# update tag 'latest'
podman tag troder/scoary-2 troder/scoary-2:latest
podman push troder/scoary-2:latest
```

## Docker / Zenodo links in Wiki

Update Zenodo:
1) Create a new release on GitHub (Title: `scoary-2:$SCOARY_VERSION`)
2) Will automatically create a new DOI on Zenodo
3) Make sure links are updated
