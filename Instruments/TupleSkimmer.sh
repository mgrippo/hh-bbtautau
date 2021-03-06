#!/bin/bash

if [ $# -ne 3 ] ; then
    echo "Usage: tree_name input_path output_path"
    exit 1
fi

TREE_NAME="$1"
INPUT_PATH="$2"
OUTPUT_PATH="$3"
INPUT_PATTERN="*.root"

mkdir -p "$OUTPUT_PATH"

IFS=$'\n'
INPUT_FILES=( $(find "$INPUT_PATH" -maxdepth 1 -type f -name "$INPUT_PATTERN" | sort) )

for input_file in "${INPUT_FILES[@]}" ; do
    input_file_name=${input_file##*/}
    output_file="$OUTPUT_PATH/$input_file_name"
    if [ -f "$output_file" ] ; then
        echo "'$input_file' was already processed."
        continue
    fi
    echo "Skimming '$input_file' -> '$output_file'..."
    ./run.sh TupleSkimmer "$TREE_NAME" "$input_file" "$output_file"
    RESULT=$?
    if [ $RESULT -ne 0 ] ; then
        echo "Error occurred while processing '$input_file'."
        exit 1
    fi
done
