#!/bin/bash

# Script to rename files from 000-1500 to have 5-digit padding
# Usage: ./rename_files.sh [directory] [file_prefix] [file_extension]
# Example: ./rename_files.sh /path/to/files cylinder .vtk

# Set default values
DIRECTORY="${1:-.}"  # Default to current directory if not specified
PREFIX="${2:-}"      # Default to empty prefix if not specified
EXTENSION="${3:-}"   # Default to empty extension if not specified

# Function to display usage
usage() {
    echo "Usage: $0 [directory] [file_prefix] [file_extension]"
    echo "Example: $0 /path/to/files cylinder .vtk"
    echo "         $0 . myfile .txt"
    echo "         $0 /some/path \"\" .dat  # for files like 001.dat, 002.dat, etc."
    echo ""
    echo "Arguments:"
    echo "  directory     - Directory containing the files (default: current directory)"
    echo "  file_prefix   - Prefix of the files to rename (can be empty)"
    echo "  file_extension - Extension of the files to rename (can be empty)"
    exit 1
}

# Check if help is requested
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    usage
fi

# Validate directory exists
if [[ ! -d "$DIRECTORY" ]]; then
    echo "Error: Directory '$DIRECTORY' does not exist!"
    exit 1
fi

echo "Renaming files in directory: $DIRECTORY"
echo "File prefix: '$PREFIX'"
echo "File extension: '$EXTENSION'"
echo ""

# Counter for renamed files
renamed_count=0
skipped_count=0

# Find the maximum number in existing files to determine the range
max_num=0
for file in "$DIRECTORY"/*; do
    if [[ -f "$file" ]]; then
        filename=$(basename "$file")
        # Extract number from filename (remove prefix and extension)
        temp_name="$filename"
        if [[ -n "$PREFIX" ]]; then
            temp_name="${temp_name#$PREFIX}"
        fi
        if [[ -n "$EXTENSION" ]]; then
            temp_name="${temp_name%$EXTENSION}"
        fi
        # Check if the remaining part is a number
        if [[ "$temp_name" =~ ^[0-9]+$ ]]; then
            num=$((10#$temp_name))  # Force base 10 to handle leading zeros
            if [[ $num -gt $max_num ]]; then
                max_num=$num
            fi
        fi
    fi
done

echo "Found files numbered up to: $max_num"
echo ""

# Loop through numbers 0 to max_num found
for i in $(seq 0 $max_num); do
    # Create the current filename patterns to search for
    # Handle different padding scenarios (1, 2, 3, 4 digits)
    for padding in {1..4}; do
        case $padding in
            1) current_num=$(printf "%d" $i) ;;
            2) current_num=$(printf "%02d" $i) ;;
            3) current_num=$(printf "%03d" $i) ;;
            4) current_num=$(printf "%04d" $i) ;;
        esac
        
        # Construct the current filename
        current_file="${PREFIX}${current_num}${EXTENSION}"
        current_path="$DIRECTORY/$current_file"
        
        # Check if this file exists
        if [[ -f "$current_path" ]]; then
            # Create the new filename with 5-digit padding
            new_num=$(printf "%05d" $i)
            new_file="${PREFIX}${new_num}${EXTENSION}"
            new_path="$DIRECTORY/$new_file"
            
            # Check if the target filename already exists and is different
            if [[ "$current_file" != "$new_file" ]]; then
                if [[ -f "$new_path" ]]; then
                    echo "Warning: Target file '$new_file' already exists. Skipping '$current_file'"
                    ((skipped_count++))
                else
                    # Perform the rename
                    mv "$current_path" "$new_path"
                    if [[ $? -eq 0 ]]; then
                        echo "Renamed: $current_file → $new_file"
                        ((renamed_count++))
                    else
                        echo "Error: Failed to rename $current_file"
                    fi
                fi
            else
                echo "File '$current_file' already has correct padding. Skipping."
                ((skipped_count++))
            fi
            break  # Found the file, no need to check other padding formats
        fi
    done
done

echo ""
echo "Renaming complete!"
echo "Files renamed: $renamed_count"
echo "Files skipped: $skipped_count"

# Make the script executable