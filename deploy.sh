#!/bin/bash

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
docker cp `dl`:/home/rstudio/_book tmp
cp -r tmp/* docs

# push changes to the website
git add docs/*
git commit -m "update the course website"
