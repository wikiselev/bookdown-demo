Channel
    .fromPath("$baseDir/course_files", type: 'dir')
    .set { ch_course_files }


Channel
    .fromPath('s3://scrnaseq-course/data/', checkIfExists: false)
    .set { ch_data }

Channel
    .fromPath('s3://scrnaseq-course/_bookdown_files/', checkIfExists: false)
    .set { ch_cached_files }

process html {

  echo true

  input: 
    file 'course_dir_work/data' from ch_data
    file 'course_dir' from ch_course_files
    path '_bookdown_files' from ch_cached_files
  
  output:
    file 'course_dir_work/website'
  
  shell:
  '''
  cp -r course_dir/* course_dir_work
  cd course_dir_work
  ln -s ../_bookdown_files .
  Rscript -e "bookdown::render_book('index.html', 'bookdown::gitbook')"
  '''
}

/*

process latex {
  
  echo true

  input:
    file 'course_dir_work/data' from ch_data2
    file 'course_dir' from ch_course_files2
    path '_bookdown_files' from ch_cached_files
  shell::
  '''
  cp -r course_dir/* course_dir_work
  cd course_dir_work
  ln -s ../_bookdown_files .
  Rscript -e "bookdown::render_book('index.html', 'bookdown::gitbook')"
  ```
}

*/
