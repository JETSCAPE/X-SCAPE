#!/usr/bin/env bash

###############################################################################
# Copyright (c) The JETSCAPE Collaboration, 2018
#
# For the list of contributors see AUTHORS.
#
# Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
#
# or via email to bugs.jetscape@gmail.com
#
# Distributed under the GNU General Public License 3.0 (GPLv3 or later).
# See COPYING for details.
###############################################################################

# Clone js-contrib into external_packages/js-contrib
# Use with: cmake -DUSE_JS_CONTRIB=ON [-DUSE_JS_FNO_HYDRO=ON] [-DUSE_JS_PYJETSCAPE=ON]

folderName="js-contrib"

if [ -d "$folderName" ]; then
  echo "$folderName already exists — skipping clone."
  exit 0
fi

git clone https://github.com/jhputschke/js-contrib.git "$folderName"
