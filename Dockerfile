FROM hemberglab/scrna.seq.course-docker

# add our scripts
ADD . /

# run scripts
CMD bash _build.sh