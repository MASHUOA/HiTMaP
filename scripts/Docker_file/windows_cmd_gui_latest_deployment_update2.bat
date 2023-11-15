cd %~dp0
docker build -t mashuoa/hitmap:gui_latest --no-cache --file Dockerfile_shiny_hitmap_upd2 .
pause