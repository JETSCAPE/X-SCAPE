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
# (MUSIC_ABORT_ON_NONFINITE=1 stops instead), see MUSIC4GPU docs/VacReset_BUG.md (PR #13),
# the in-memory evolution store releasing its memory when an event is cleaned
# instead of keeping its largest size for the whole run (MUSIC4GPU PR #15),
# and less host CPU per event (MUSIC4GPU PR #18): MUSIC_CUDA_SYNC=block|yield|spin|auto
# selects how the host waits for the GPU (unset keeps CUDA's default), and the string
# source computes per-string constants once per step and skips zero envelopes and,
# without preflow, the transverse-flow terms (output bit-identical),
# and less host work per run (MUSIC4GPU PR #19): no source fills for charges that
# are not evolved, the output_momentum_anisotropy_vs_etas switch, frames packed in
# place into a huge-page store (pinned staging on discrete GPUs), the no-surface
# edge test as a GPU reduction, and store_hydro_info_frames (output bit-identical),
# and MusicCheckpoint (MUSIC4GPU music2_restart): a run starts from another run's
# saved state, the jet leg of the two-stage chain from the background's
commitHash="13030792213a9d1eec092402f1a8c93e0bc4154f"

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