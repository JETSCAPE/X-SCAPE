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

# clone the lbt-tables repository
repo_dir="LBT-tables-repo"
git clone --depth 1 https://github.com/JETSCAPE/LBT-tables.git "$repo_dir"

# temporary directory where the split LBT tables will be extracted
split_extract_dir="tmp/LBT-tables-split-extract"

# where the final unsplit LBT tables will be placed
unsplit_dir="LBT-tables"

if [ ! -d "$repo_dir/LBT-tables-split-archive" ]; then
    echo "Error archive directory does not exist."
    exit 1
fi

mkdir -p $split_extract_dir
mkdir -p $unsplit_dir

# extract all the .tar.gz files to the temporary directory
for file in "$repo_dir/LBT-tables-split-archive"/*.tar.gz; do
    if [ -f "$file" ]; then
        tar -xzf "$file" -C "$split_extract_dir"
        echo "Extracted '$file' to '$split_extract_dir'"
    else
        echo "Error '$file' does not conform to .tar.gz format."
        exit 1
    fi
done

# the unique base names: myFile.part0 and myFile.part1 share the base name myFile.
# A sorted list instead of an associative array, which the bash 3.2 of macOS lacks.
base_names=$(
    for file in "$split_extract_dir"/*; do
        filename=$(basename "$file")
        if [[ "$filename" =~ \.part([0-9]+)$ ]]; then
            echo "${filename%.part*}"
        else
            echo "Error: '$filename' - all files must end with .partN, where N is a positive integer." >&2
        fi
    done | sort -u
)

while IFS= read -r base_name; do
    [ -n "$base_name" ] || continue
    output_file="$unsplit_dir/$base_name"

    # if a file is small enough that it has only one part (.part0), copy it to the unsplit directory
    if [ -f "$split_extract_dir/$base_name.part0" ] && [ ! -f "$split_extract_dir/$base_name.part1" ]; then
        cp "$split_extract_dir/$base_name.part0" "$output_file"
        echo "Copied '$split_extract_dir/$base_name.part0' to '$output_file'"
    else
        # concatenate the parts in numerical order (.part2 before .part10), counting up from
        # .part0; not with sort -V, which orders them as text on macOS
        > "$output_file"
        i=0
        while [ -f "$split_extract_dir/$base_name.part$i" ]; do
            cat "$split_extract_dir/$base_name.part$i" >> "$output_file"
            echo "Added '$split_extract_dir/$base_name.part$i' to '$output_file'"
            i=$((i + 1))
        done
    fi
done <<< "$base_names"

# remove the temporary extracted split files
rm -rf tmp

# remove the local lbt-tables repository if not wanted
# rm -rf "$repo_dir"
