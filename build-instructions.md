## Instructions for Building Course 

### Clone Repository
In order to download the course. Enter the following command into the directory you want the course to be downloaded into:
```
git clone https://github.com/cellgeni/scrnaseq-course-private.git
```

### Installing the Image
The course uses a docker image within a singularity environment. In order to ensure you have all the correct
dependencies installed please download the v4.07 [docker image](https://quay.io/repository/hemberg-group/scrna-seq-course?tab=tags).

A specific [version of singularity](https://github.com/hpcng/singularity/tree/v3.5.3) is needed.

There are also [instructions](https://github.com/hpcng/singularity/blob/v3.5.3/INSTALL.md) for installing singularity.

The software nextflow is also used to build the course which has its own [installation instructions](https://www.nextflow.io/docs/latest/getstarted.html).

### How to Build the Course
In order to build the course and generate new cache files please input the following code into a file (i.e. run-course):
```
vi run-course
```
then copy the following code into the file:
```
#!/bin/bash

source=cellgeni/scrnaseq-course-private
source=/path-to-current-directory/scrnaseq-course-private/main.nf

export PATH=$PATH:/path-to-installed-singularity-software/singularity-v3.5.3/bin/

set -euo pipefail

nextflow run $source -profile singularity -with-report reports/report.html -resume -ansi-log false
```

To then run this code, use the following command:
```
/path-to-directory-containing-run-course-file/run-course
```

Or if the file is in your current directory then you can use:
```
./run-course
```

This should build the course. The work directory will be provided at the end and the newly built
cache files will be located at:
```
/path-to-work-dir/course_work_dir/_bookdown_files/
```

### How to Upload Newly Built Cache to Amazon S3 Bucket
If you need to upload new cache to the github repo then you will need the AWS Access Key ID and 
AWS Secret Access Key (not provided here).

Then you need to start a singularity shell using the following command:
```
SINGULARITYENV_AWS_ACCESS_KEY_ID=NOT-PROVIDED  \
SINGULARITYENV_AWS_SECRET_ACCESS_KEY=STILL-NOT-PROVIDED  \
/path-to-installed-singularity-software/singularity-v3.5.3/bin/singularity shell -B /any-paths-that-need-to-be mounted /path-to-docker-images/quay.io-hemberg-group-scrna-seq-course-v4.07.img
```

Once you have the shell started use the following command to upload new cache:
```
aws s3 sync /path-to-work-dir/course_work_dir/_bookdown_files/ s3://scrnaseq-course/_bookdown_files/
```
