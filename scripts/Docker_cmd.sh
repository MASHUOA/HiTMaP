docker pull mashuoa/hitmap
docker tag mashuoa/hitmap:1.02 mashuoa/hitmap:latest

docker ps
docker commit cf8d5717bbbf mashuoa/hitmap:1.02
docker image ls
docker image push mashuoa/hitmap:1.02