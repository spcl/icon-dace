import yaml
import sys
import subprocess

import dace

def read_yaml(file_path):
    with open(file_path, 'r') as file:
        return yaml.safe_load(file)

def process_entries(entries, sdfg_sources_path, src_dir):
    for file_path, sdfg_names in entries.items():
        print(f"File: {file_path}")
        for sdfg_name, lines in sdfg_names.items():
            print(f"  SDFG Name: {sdfg_name}")
            command = [
                'python', f'{sdfg_sources_path}/generate.py',
                f'{sdfg_name}',
                f'{sdfg_sources_path}/integrations.yaml',
                f'{sdfg_sources_path}'
            ]
            print(f"Generating SDFG {sdfg_name}")
            print(f"Executing command: {' '.join(command)}")
            subprocess.run(command, check=True)

            command = [
                'python', f'{sdfg_sources_path}/optimize_{sdfg_name}.py',
                f'{sdfg_sources_path}/{sdfg_name}_unsimplified.sdfgz',
                f'{src_dir}/{sdfg_name}_optimized.sdfgz'
            ]
            print(f"Optimizing {sdfg_name}, writing output to: {src_dir}/{sdfg_name}_optimized.sdfgz")
            print(f"Executing command: {' '.join(command)}")
            subprocess.run(command, check=True)

            command = [
                'python', f'{sdfg_sources_path}/optimize_{sdfg_name}.py',
                f'{sdfg_sources_path}/{sdfg_name}_unsimplified.sdfgz',
                f'{src_dir}/{sdfg_name}_optimized.sdfgz'
            ]
            print(f"Compiling {src_dir}/{sdfg_name}_optimized.sdfgz")
            print(f"Executing command: {' '.join(command)}")
            subprocess.run(command, check=True)

            sdfg = dace.SDFG.from_file(f'{src_dir}/{sdfg_name}_optimized.sdfgz')
            sdfg.compile()

            for line_range in lines:
                print(f"    Lines: {line_range}")
                # Call further Python scripts based on the input
                # Example: subprocess.call(['python', 'your_script.py', file_path, function, str(line_range)])

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <src_dir> <path_to_icon_sources> <path_to_integrations.yaml")
        sys.exit(1)

    src_dir = sys.argv[1]
    sdfg_sources_path = sys.argv[2]
    yaml_file_path = sys.argv[3]
    entries = read_yaml(yaml_file_path)
    process_entries(entries, sdfg_sources_path, src_dir)