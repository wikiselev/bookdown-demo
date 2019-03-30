#!/usr/bin/env bash

# download course data
aws --endpoint-url https://cog.sanger.ac.uk --no-sign-request s3 cp s3://scrnaseq-course/data data --recursive