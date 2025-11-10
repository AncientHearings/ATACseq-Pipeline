#!/bin/bash

success=true

for dir in results/*; do
    if [ -d "$dir" ]; then
        echo "......................."
        echo "Directory: ${dir}"
        echo "Creating .gitkeep file"
        if touch "${dir}/.gitkeep"; then
            echo "File created"
        else
            echo "ERROR: Failed to create file in ${dir}"
            success=false
        fi
        echo "......................."
    fi
done

if [ "$success" = true ]; then
    echo "Successful execution of the for loop"
else
    echo "For loop encountered errors."
fi

echo "Done"
