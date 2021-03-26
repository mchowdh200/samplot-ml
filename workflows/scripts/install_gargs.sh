#!/bin/env bash
gargs=$1
install_dir=$(dirname $gargs)
[[ ! -d $install_dir ]] && mkdir $install_dir
case "$(uname -s)" in
    Darwin)
        wget https://github.com/brentp/gargs/releases/download/v0.3.9/gargs_darwin \
             -O $gargs
        ;;
    Linux)
        wget https://github.com/brentp/gargs/releases/download/v0.3.9/gargs_linux \
             -O $gargs
        ;;
esac
chmod +x $gargs
