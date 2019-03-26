
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
    .fromPath('*')
    .set { ch_course_files }

process foo {
  input: 
    file fs from ch_course_files.collect()
  script:
  """
  Rscript -e "bookdown::render_book('index.html', 'bookdown::gitbook')"
  """
}