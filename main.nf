
// process foo {
//   input:
//   val x from 1
//   output:
//   file 'x.txt' into result

//   """
//   echo $x > x.txt
//   """
// }

Channel
    .fromPath('course_files', type: 'dir')
    .into { ch_course_files1, ch_course_files2 }

process html {
  input: 
    file fs from ch_course_files1
  script:
  """
  ln -s course_files/* .
  Rscript -e "bookdown::render_book('index.html', 'bookdown::gitbook')"
  """
}

process latex {
  input: 
    file fs from ch_course_files2
  script:
  """
  ln -s course_files/* .
  Rscript -e "bookdown::render_book('index.html', 'bookdown::pdf_book')"
  """
}
