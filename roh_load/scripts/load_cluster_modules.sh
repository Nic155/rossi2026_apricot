#!/usr/bin/env bash

# Source this file from submit_pipeline.sh and the SLURM job scripts.
# It loads the cluster modules required by this pipeline.

if ! command -v module >/dev/null 2>&1; then
    for module_init in \
        /etc/profile.d/modules.sh \
        /usr/share/Modules/init/bash \
        /usr/share/lmod/lmod/init/bash
    do
        if [[ -r "$module_init" ]]; then
            # shellcheck disable=SC1090
            source "$module_init"
            break
        fi
    done
fi

if ! command -v module >/dev/null 2>&1; then
    printf 'ERROR: Environment Modules command "module" is not available.\n' >&2
    return 1 2>/dev/null || exit 1
fi

module load devel/python/Python-3.11.1
module load bioinfo/bedtools/2.31.1
