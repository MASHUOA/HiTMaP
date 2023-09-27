cd %~dp0
docker build -t mashuoa/hitmap:gui_latest --file Dockerfile_shiny_hitmap .
pause