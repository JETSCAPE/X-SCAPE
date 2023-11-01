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
##############################################################################

# download the code package
commitHash="4bc2badcd31401bcdea8be5a2efc778bdf99fc57" # for xscape 1.1
rm -fr iSS
git clone https://github.com/chunshen1987/iSS -b JETSCAPE iSS
cd iSS
git checkout $commitHash
#git checkout tags/v1.1.1 -b v1.1.1
rm -fr iSS/.git
