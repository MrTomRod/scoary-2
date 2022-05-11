# How to publish

## Pypi

Create update:

1) Change version number in [pyproject.toml](pyproject.toml)
2) Create new package and upload:

```bash
# build: will create files in dist/
poetry build
# test: install .whl file
pip install dist/scoary_2-*-py3-none-any.whl
# upload
poetry publish
```

## Docker / Podman

If you use docker, simply replace each `podman` with `docker`.

```shell
podman build --tag troder/scoary-2 .
```

Publish docker image:

```shell
# podman login docker.io --get-login
# podman login docker.io
podman tag troder/scoary-2 troder/scoary-2:<scoary-version>
podman push troder/scoary-2:<scoary-version>

# update tag 'latest'
podman tag troder/scoary-2 troder/scoary-2:latest
podman push troder/scoary-2:latest
```
