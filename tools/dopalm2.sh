#!/bin/bash

set -eu

TMPDIR="$(mktemp -d)"
trap 'rm -r -- "${TMPDIR}"' EXIT

palm2() {

    while getopts ":f:" opt; do
        case $opt in
            f)
                local flag=$OPTARG
                ;;
            \?)
                echo "Invalid option: -$OPTARG" >&2
                ;;
            :)
                echo "Option -$OPTARG requires an argument." >&2
                exit 1
        esac
    done
    shift "$((OPTIND -1))"

    # gather files
    local niis=("$@")
    
    local merged_series
    merged_series=$(mktemp --suffix=.nii.gz -p "${TMPDIR}")
    local mask=/fastscratch/myscratch/pssadil/meta/tools/MNI152_T1_2mm_brain_mask.nii
    
    fslmerge -t ${merged_series} "${niis[@]}"
    gunzip "${merged_series}"
    infile="${merged_series%.gz}"
    
    palm -i ${infile} -m ${mask} -o "${flag}" -T -logp -ise -saveglm -savedof -savemetrics -saveparametric -savemax

}
export -f palm2

palm2 "$@"
