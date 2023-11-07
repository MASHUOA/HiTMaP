cd %~dp0
docker build -t mashuoa/hitmap:latest --file Dockerfile_base_latest_upd . 2> log.txt
pause