#!/bin/bash

md=$1
html=$2
header=$3
nav=$4
pandoc_toc_sidebar=$5


cd \
${pandoc_toc_sidebar} \
&& \
pandoc \
--self-contained \
--toc \
--highlight-style=haddock \
-f markdown \
-t html5 \
--include-in-header ${header} \
--template=toc-sidebarL.html \
-B ${nav} \
${md} \
> ${html}