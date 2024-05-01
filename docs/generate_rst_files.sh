#!/bin/bash

# Directory containing your .py files
SOURCE_DIR="../src/bakir"
# Directory to store the generated .rst files
DOCS_DIR="./"

# Check if the documentation directory exists, create if not
if [ ! -d "$DOCS_DIR" ]; then
    mkdir -p "$DOCS_DIR"
fi

# Loop through all .py files in the source directory
for file in $SOURCE_DIR/*.py; do

    if [[ "$file" == *"data"* ]]; then
        echo "Skipping $file as it contains 'data'"
        continue

    # Extract the filename without the extension
    filename=$(basename -- "$file")
    module_name="${filename%.*}"

    # Skip __init__.py files
    if [ "$module_name" == "__init__" ]; then
        continue
    fi

    # Create an .rst file with the same name in the documentation directory
    rst_file="$DOCS_DIR/$module_name.rst"
    echo "Creating $rst_file..."

    # Write contents to the .rst file
    echo "$module_name" > "$rst_file"
    echo "=====================" >> "$rst_file"
    echo "" >> "$rst_file"
    echo ".. automodule:: bakir.$module_name" >> "$rst_file"
    echo "   :members:" >> "$rst_file"
    echo "   :undoc-members:" >> "$rst_file"
    echo "   :show-inheritance:" >> "$rst_file"
done

echo "Documentation files have been generated."
