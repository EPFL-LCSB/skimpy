#!/bin/bash

#wget https://dl.google.com/linux/direct/google-chrome-stable_current_amd64.deb


function install_chrome_browser() {
    echo '>>> Installing Chrome'

    wget -q -O - https://dl-ssl.google.com/linux/linux_signing_key.pub | apt-key add -

    echo "deb [arch=amd64] http://dl.google.com/linux/chrome/deb/ stable main" >> /etc/apt/sources.list.d/google.list

    apt-get update
    apt-get install -y google-chrome-stable --no-install-recommends

    # Disable sandboxing - it conflicts with unprivileged lxc containers
    sed -i 's|HERE/chrome"|HERE/chrome" --disable-setuid-sandbox --enable-logging --no-sandbox|g' \
               "/opt/google/chrome/google-chrome"
}

install_chrome_browser

#apt install -y ./google-chrome-stable_current_amd64.deb

