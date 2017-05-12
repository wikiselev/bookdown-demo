#!/bin/bash

WORKSPACE=$1

# Delete every Docker containers
# Must be run first because images are attached to containers
docker rm -f $(docker ps -a -q)

# Delete every Docker image
docker rmi -f $(docker images -q)

# get the docker
docker pull quay.io/hemberg-group/scrna-seq-course:latest
# run the docker
docker run quay.io/hemberg-group/scrna-seq-course

# copy files from the docker
alias dl='docker ps -l -q'
docker cp `dl`:_book $WORKSPACE/tmp1
cp -r tmp1/* docs
docker cp `dl`:tung $WORKSPACE/tmp2
cp -r tmp2/* tung

# push changes to the website
git add docs/*
git add tung/*
git commit -m "update the course website"
# git push origin master
