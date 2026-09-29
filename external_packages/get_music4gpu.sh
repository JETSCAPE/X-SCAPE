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

# using a commit from the MUSIC repository that is compatible with the current X-SCAPE version
folderName="music4gpu"
# XSCAPE branch with the jet source slot (add_hydro_source_terms_from_jet, which
# MusicWrapper calls since X-SCAPE PR #138), per-step jet droplet pruning, the
# freeze_out_surface switch, and the GPU-path source fill speed-ups (MUSIC4GPU PR #10:
# skip steps where no source can deposit -- HydroSourceJETSCAPE reports the jet side
# since X-SCAPE PR #144 -- and bin the strings by transverse reach; output bit-identical),
# and StringFind4 failing loudly on a parameter file without EndOfData instead of
# hanging (MUSIC4GPU PR #11), the parallel, deterministic freeze-out surface search
# (MUSIC4GPU PR #12), and the GPU fix for grids freezing in dilute regions: vacuum
# cells at rest, a guard against non-finite W^{mu nu}/Pi with a counter and warning
# (MUSIC_ABORT_ON_NONFINITE=1 stops instead), see VacReset_BUG.md (MUSIC4GPU PR #13)
commitHash="15ec5e3ba87fa89749c1f06aad479add4aa0f556"

git clone https://github.com/jhputschke/MUSIC4GPU.git -b XSCAPE $folderName
cd $folderName
git checkout $commitHash
cd EOS
bash download_hotQCD.sh SMASH_binary
# EOS 9 (UrQMD hadron list), used by e.g. the PyJetscape prod_AuAu_0_10 productions
bash download_hotQCD.sh binary

# Finite muB not supprted on GPU yet!!!
# Download the 4D EoS tables (only needed for EOS 20, only download if necessary, large files)
# bash download_Neos4D.sh UrQMD