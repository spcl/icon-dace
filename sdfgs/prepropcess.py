import sys

def replace_in_file(file_path):
    try:
        # Read the content of the file
        with open(file_path, 'r', encoding='utf-8') as file:
            content = file.read()

        # Replace occurrences
        content = content.replace("geofac_qdiv", "geofac_div")
        content = content.replace("cells_plwa_verts", "cells_aw_verts")
        content = content.replace("ptr_patch", "p_patch")

        # Write the modified content back to the file
        with open(file_path, 'w', encoding='utf-8') as file:
            file.write(content)

        print(f"Replacements made successfully in {file_path}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <file_path>")
        sys.exit(1)

    file_path = sys.argv[1]
    replace_in_file(file_path)
