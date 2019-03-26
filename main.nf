
// process foo {
//   input:
//   val x from 1
//   output:
//   file 'x.txt' into result

//   """
//   echo $x > x.txt
//   """
// }

ch_course_files = Channel.fromPath('course_files', type: 'dir')

process foo {
  input: 
    file fs from ch_course_files
  script:
  """
  ln -s course_files/* .
  Rscript -e "bookdown::render_book('index.html', 'bookdown::gitbook')"
  """
}
