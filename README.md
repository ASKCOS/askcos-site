# askcos-site
Web application for the prediction of feasible synthetic routes towards a desired compound and associated tasks related to synthesis planning. Originally developed under the DARPA Make-It program and now being developed under the [MLPDS Consortium](http://mlpds.mit.edu).

## 2020.07 Release Notes

User notes:

Developer notes:

Bug fixes:

For old release notes, see the [ASKCOS releases page](https://gitlab.com/mlpds_mit/ASKCOS/ASKCOS/-/releases).

## Getting Started

This package provides a web interface for the [`askcos-core`](https://gitlab.com/mlpds_mit/ASKCOS/askcos-core) Python package. It is built using Django, with newer pages using Vue.js for dynamic content. Celery is also used for asyncronous task processing. We recommend deployment using Docker, for which scripts can be found in the [`askcos-deploy`](https://gitlab.com/mlpds_mit/ASKCOS/askcos-deploy) repository.

### Downloading with GitLab Deploy Tokens

This repository can be downloaded using deploy tokens, which provide __read-only__ access to the source code and our container registry in GitLab. The deploy tokens can be found on the [MLPDS Member Resources ASKCOS Versions Page](https://mlpds.mit.edu/member-resources-releases-versions/). The only software prerequisites are git, docker, and docker-compose.

```bash
$ export DEPLOY_TOKEN_USERNAME=
$ export DEPLOY_TOKEN_PASSWORD=
$ git clone https://$DEPLOY_TOKEN_USERNAME:$DEPLOY_TOKEN_PASSWORD@gitlab.com/mlpds_mit/askcos/askcos-site.git
```

### Building a Docker Image

Before building the `askcos-site` image, you log in to the ASKCOS GitLab registry to download the `askcos-core` image which is used as the base image (or build it yourself). Use the same deploy tokens as above to log in to the registry:

```bash
docker login registry.gitlab.com -u $DEPLOY_TOKEN_USERNAME -p $DEPLOY_TOKEN_PASSWORD
```

Then, the askcos-site image can be built using the Dockerfile in this repository.

```bash
$ cd askcos-site
$ docker build -t <image name> .
```

A Makefile is also provided to simplify the build command by providing a default image name and tag:

```bash
$ cd askcos-site
$ make build
```