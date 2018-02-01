FROM quay.io/hemberg-group/scrna-seq-course-r:latest

# add our scripts
ADD . /home/rstudio/

# run scripts
CMD cd /home/rstudio && bash build.sh
