#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 input_folder output_folder command"
    exit 1
fi

INPUT_FOLDER="$1"
OUTPUT_FOLDER="$2"
COMMAND="$3"

# Function to process files
process_files() {
    local input_dir="$1"
    local output_dir="$2"

    # Create the output directory if it doesn't exist
    mkdir -p "$output_dir"

    # Loop through each item in the input directory
    for item in "$input_dir"/*; {
        if [ -d "$item" ]; then
            # If it's a directory, recurse into it
            local subdir=$(basename "$item")
            process_files "$item" "$output_dir/$subdir"
        elif [ -f "$item" ]; then
            # If it's a file, process it with the command line tool
            local filename=$(basename "$item")
            local output_file="$output_dir/$filename"
            $COMMAND "$item" "$output_file"
        fi
    }
}

# Ensure the output directory is clean
rm -rf "$OUTPUT_FOLDER"
mkdir -p "$OUTPUT_FOLDER"

# Start processing files from the input directory
process_files "$INPUT_FOLDER" "$OUTPUT_FOLDER"

