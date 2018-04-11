#!/usr/bin/env bash

export DYLD_LIBRARY_PATH=""
epstopdf $1 --nocompress --autorotate=All
