docker pull mashuoa/hitmap
docker tag mashuoa/hitmap:1.02 mashuoa/hitmap:latest

docker ps
docker commit cf8d5717bbbf mashuoa/hitmap:1.02
docker image ls
docker image push mashuoa/hitmap:1.02
docker build --no-cache -t mashuoa/hitmap:1.06 -f C:\Users\admgguo484\Documents\GitHub\HiTMaP\dockerfiles_shiny_hitmap
docker build --no-cache -t mashuoa/hitmap:1.06 .
docker run --rm -p 80:3838 -v %userprofile%\Documents\expdata:/root/expdata mashuoa/hitmap

ssh gguo484@md-cer00340.its.auckland.ac.nz
password:2authentoken
