pipeline {
  agent any
  stages {
    stage('Clean up old stuff') {
      steps {
        sh '''# Delete old Docker containers and images
docker rm -f $(docker ps -a -q)
docker rmi -f $(docker images -q)'''
      }
    }
    stage('Run Docker') {
      steps {
        sh 'docker run quay.io/hemberg-group/scrna-seq-course:latest'
      }
    }
    stage('Copy from Docker') {
      steps {
        sh '''# copy files from the docker
alias dl=\'docker ps -l -q\'
docker cp `dl`:/home/rstudio/_book $WORKSPACE/tmp1
cp -r tmp1/* docs'''
      }
    }
    stage('Commit changes') {
      steps {
        sh '''# commit changes
git add docs/*
git commit -m "update the course website"
git push origin master
'''
      }
    }
  }
  triggers {
    cron('H 14 * * *')
  }
}