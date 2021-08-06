# askcos-site
[![askcos-base](https://img.shields.io/badge/-askcos--base-lightgray?style=flat-square)](https://github.com/ASKCOS/askcos-base)
[![askcos-data](https://img.shields.io/badge/-askcos--data-lightgray?style=flat-square)](https://github.com/ASKCOS/askcos-data)
[![askcos-core](https://img.shields.io/badge/-askcos--core-lightgray?style=flat-square)](https://github.com/ASKCOS/askcos-core)
[![askcos-site](https://img.shields.io/badge/-askcos--site-blue?style=flat-square)](https://github.com/ASKCOS/askcos-site)
[![askcos-deploy](https://img.shields.io/badge/-askcos--deploy-lightgray?style=flat-square)](https://github.com/ASKCOS/askcos-deploy)

Web application for the prediction of feasible synthetic routes towards a desired compound and associated tasks related to synthesis planning. Originally developed under the DARPA Make-It program and now being developed under the [MLPDS Consortium](http://mlpds.mit.edu).

## Getting Started

This package provides a web interface for the [`askcos-core`](https://github.com/ASKCOS/askcos-core) Python package. It is built using Django, with newer pages using Vue.js for dynamic content. Celery is also used for asyncronous task processing. We recommend deployment using Docker, for which scripts can be found in the [`askcos-deploy`](https://github.com/ASKCOS/askcos-deploy) repository. Additional information on deployment and releases can be found at [askcos.gitlab.io](https://askcos.gitlab.io).

### Building a Docker Image

The `askcos-site` image can be built using the Dockerfile in this repository. It depends on Docker images for `askcos-data` and `askcos-site`, which can be built manually or pulled from Docker Hub.

```bash
$ cd askcos-site
$ docker build -t <image name> .
```

A Makefile is also provided to simplify the build command by providing a default image name and tag:

```bash
$ cd askcos-site
$ make build
```
