#!/bin/bash

docker build -t hscc25 .
docker run --memory=16g --rm -v "${PWD}/output:/app/output" hscc25