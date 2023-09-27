cd %~dp0
docker build -t mashuoa/hitmap:latest --no-cache --file Dockerfile_base_latest .
pause