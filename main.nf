
Channel
    .fromPath('course_files', type: 'dir')
    .into { ch_course_files1; ch_course_files2 }

Channel
    .fromPath('s3://scrnaseq-course.cog.sanger.ac.uk/data', type: 'dir')
    .into { ch_data1; ch_data2 }

process html {
  input: 
    file fs from ch_course_files1
    file dat from ch_data1
  script:
  """
  cp -r course_files/* .
  Rscript -e "bookdown::render_book('index.html', 'bookdown::gitbook')"
  """
}

process latex {
  input: 
    file fs from ch_course_files2
  script:
  """
  cp -r course_files/* .
  Rscript -e "bookdown::render_book('index.html', 'bookdown::pdf_book')"
  """
}
