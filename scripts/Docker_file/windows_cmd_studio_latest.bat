cd %~dp0
docker build -t mashuoa/hitmap:studio --no-cache --file Dockerfile_studio_latest2 . 2> log.txt
pause