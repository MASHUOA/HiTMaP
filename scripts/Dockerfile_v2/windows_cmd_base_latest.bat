cd %~dp0
docker build --no-cache -t mashuoa/hitmap:base_prep --file Dockerfile_base_prep . 2> log_base_prep.txt
docker build --no-cache -t mashuoa/hitmap:latest --file Dockerfile_base_latest . 2> log_base.txt
pause