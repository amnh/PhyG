#!/bin/bash
set -e

TUTORIAL_STEM='tutorial'
TUTORIAL_DIR="${TUTORIAL_STEM}s"
TUTORIAL_DOCS_PATTERN="${TUTORIAL_STEM}_*_doc"

note() {
    local above="┍━━"
    local below="┕━━"
    local prefix="│ "
    local pattern="${above}\n%s\n${below}\n"
    echo -n -e "${1}" | sed "s/^/$prefix /" | xargs -0 printf ${pattern}
}

build_tutorial_documetation() {
    local DOC_PATH="${1}"
    local DIR_INIT=$(pwd)
    cd "${DOC_PATH}"

    # Get file LaTeX file
    FILE_TEX=$(find . -maxdepth 1 -iname "*.tex" -print -quit)
    FILE_BASENAME="${FILE_TEX%.*}"

    # Build PDF from LaTeX
    latexmk -f -pdf "${FILE_BASENAME}" > /dev/null 2>&1
    ret=$?
    if [ $ret -ne 0 ]; then
        printf "\n>>>\tFailed to build PDF for '%s'!\n>>>\tCheck the log file:>>>\t\t%s.log" "${FILE_BASENAME}" "${FILE_BASENAME}"
    fi

    # Clean up build artifacts and exit
    latexmk -c "${FILE_BASENAME}" > /dev/null 2>&1
    cd "${DIR_INIT}"
}

# Find all tutorials
TUTORIAL_FOUND=$(find ${TUTORIAL_DIR} -maxdepth 2 -iname "${TUTORIAL_DOCS_PATTERN}" -print | sort)

# Conert found results to a Bash "array" type
SAVEIFS=${IFS} # Save current IFS (Internal Field Separator)
IFS=$'\n'      # Change IFS to newline char
TUTORIAL_FILEPATH=(${TUTORIAL_FOUND})
IFS=${SAVEIFS} # Restore original IFS


printf "Compiling LaTeX documentation for the '${#TUTORIAL_FILEPATH[@]}' tutorials found:\n\n"
for (( i=0; i<${#TUTORIAL_FILEPATH[@]}; i++ ))
do
    printf "\t%d:  %s\n" "${i}" "${TUTORIAL_FILEPATH[${i}]}"
    build_tutorial_documetation "${TUTORIAL_FILEPATH[${i}]}"
done
printf "Documentation PDFs successfully built!\n
