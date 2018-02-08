pipeline {
  agent any
  stages {
    stage('deploy') {
      steps {
        sh '''# build and deploy
sh deploy.sh $WORKSPACE
'''
      }
    }
  }
}