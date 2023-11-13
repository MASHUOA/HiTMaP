cd %~dp0
docker build --no-cache -t mashuoa/hitmap:latest --file Dockerfile_base_latest . 2> log_base.txt
pause