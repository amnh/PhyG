#!/bin/bash
set -e

TUTORIAL_STEM='PhyG-User-Manual'
TUTORIAL_DIR="doc/${TUTORIAL_STEM}"

note() {
    local above="┍━━"
    local below="┕━━"
    local prefix="│ "
    local pattern="${above}\n%s\n${below}\n"
    echo -n -e "${1}" | sed "s/^/$prefix /" | xargs -0 printf ${pattern}
}

DIR_INIT=$(pwd)
cd "${TUTORIAL_DIR}"

# Get file LaTeX file
FILE_TEX=$(find . -maxdepth 1 -iname "*.tex" -print -quit)
FILE_BASENAME="${FILE_TEX%.*}"

printf "FILE_BASENAME:\t%s\n" "${FILE_BASENAME}"

# Build PDF from LaTeX
latexmk -f -pdf "${FILE_BASENAME}" > /dev/null 2>&1
ret=$?
if [ $ret -ne 0 ]; then
    printf "\n>>>\tFailed to build PDF for '%s'!\n>>>\tCheck the log file:>>>\t\t%s.log" "${FILE_BASENAME}" "${FILE_BASENAME}"
    exit 1
fi

# Clean up build artifacts and exit
latexmk -c "${FILE_BASENAME}" > /dev/null 2>&1
cd "${DIR_INIT}"

printf "User Manual PDF successfully built!\n
