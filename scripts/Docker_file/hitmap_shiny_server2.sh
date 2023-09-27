# Set up directories and permissions
if [ -x "$(command -v rstudio-server)" ]; then
    DEFAULT_USER=${DEFAULT_USER:-rstudio}
    adduser "${DEFAULT_USER}" shiny
fi

cp -R /usr/local/lib/R/site-library/HiTMaP/HiTMaP_GUI/* /srv/shiny-server
chown shiny:shiny /var/lib/shiny-server
mkdir -p /var/log/shiny-server
chown shiny:shiny /var/log/shiny-server

# create init scripts
mkdir -p /etc/services.d/shiny-server
cat <<"EOF" >/etc/services.d/shiny-server/run
#!/usr/bin/with-contenv bash
## load /etc/environment vars first:
for line in $( cat /etc/environment ) ; do export $line > /dev/null; done
if [ "$APPLICATION_LOGS_TO_STDOUT" != "false" ]; then
    exec xtail /var/log/shiny-server/ &
fi
exec shiny-server 2>&1
EOF
chmod +x /etc/services.d/shiny-server/run

# install init script
cp /rocker_scripts/init_set_env.sh /etc/cont-init.d/01_set_env

# Clean up
rm -rf /var/lib/apt/lists/*
rm -rf /tmp/downloaded_packages

## Strip binary installed lybraries from RSPM
## https://github.com/rocker-org/rocker-versioned2/issues/340
chmod 777 /usr/local/lib/R/site-library/*
strip /usr/local/lib/R/site-library/*/libs/*.so
