s3cmd -r put website s3://scrnaseq-course/
s3cmd setacl --acl-public s3://scrnaseq-course
s3cmd setacl --acl-public -r s3://scrnaseq-course
s3cmd modify --mime-type text/html s3://scrnaseq-course/website/*.html
s3cmd modify --mime-type text/css s3://scrnaseq-course/website/*.css
s3cmd modify --mime-type text/css s3://scrnaseq-course/website/libs/gitbook-2.6.7/css/*.css
