import yaml
import sys
import subprocess
import os
import dace


def read_yaml(file_path):
    with open(file_path, "r") as file:
        return yaml.safe_load(file)


def get_folder_and_filename(file_path):
    folder_path = os.path.dirname(file_path)
    filename = os.path.basename(file_path)
    filename_without_ext = os.path.splitext(filename)[0]
    return folder_path, filename_without_ext


def process_entries(entries, sdfg_sources_path, src_dir):
    for file_path, sdfg_names in entries.items():
        print(f"File: {file_path}")
        for sdfg_name, lines in sdfg_names.items():
            print(f"  SDFG Name: {sdfg_name}")
            command = [
                "python",
                f"{sdfg_sources_path}/generate.py",
                f"{sdfg_name}",
                f"{sdfg_sources_path}/integrations.yaml",
                f"{sdfg_sources_path}",
            ]
            print(f"Generating SDFG {sdfg_name}")
            print(f"Executing command: {' '.join(command)}")
            subprocess.run(command, check=True)

            command = [
                "python",
                f"{sdfg_sources_path}/optimize_{sdfg_name}.py",
                f"{sdfg_sources_path}/{sdfg_name}_unsimplified.sdfgz",
                f"{sdfg_sources_path}/{sdfg_name}_optimized.sdfgz",
            ]
            print(
                f"Optimizing {sdfg_name}, writing output to: {sdfg_sources_path}/{sdfg_name}_optimized.sdfgz"
            )
            print(f"Executing command: {' '.join(command)}")
            subprocess.run(command, check=True)

            command = [
                "python",
                f"{sdfg_sources_path}/optimize_{sdfg_name}.py",
                f"{sdfg_sources_path}/{sdfg_name}_unsimplified.sdfgz",
                f"{sdfg_sources_path}/{sdfg_name}_optimized.sdfgz",
            ]
            print(f"Compiling {sdfg_sources_path}/{sdfg_name}_optimized.sdfgz")
            print(f"Executing command: {' '.join(command)}")
            subprocess.run(command, check=True)

            sdfg = dace.SDFG.from_file(
                f"{sdfg_sources_path}/{sdfg_name}_optimized.sdfgz"
            )
            sdfg.compile()

            module_filepath = file_path
            folder_name, module_name = get_folder_and_filename(module_filepath)
            assert module_filepath.endswith(".f90") or module_filepath.endswith(
                ".F90"
            ), f"{module_filepath}"
            suffix = module_filepath.split(".")[-1]

            command = [
                "python",
                f"{sdfg_sources_path}/utils/dace-genfi.py",
                f"{sdfg_sources_path}/{sdfg_name}_optimized.sdfgz",
                f"{sdfg_sources_path}/meta_data.yaml",
                f"{sdfg_sources_path}/{sdfg_name}_module_definitions.yaml",
                f"{src_dir}/{folder_name}/{sdfg_name}_bindings.{suffix}",
            ]
            print(f"Executing command: {' '.join(command)}")
            subprocess.run(command, check=True)

            command = [
                "python",
                f"{sdfg_sources_path}/utils/dace-substf.py",
                f"{sdfg_sources_path}/integrations.yaml",
                f"{sdfg_sources_path}/{sdfg_name}_associations.yaml",
                f"{src_dir}/{file_path}",
                f"{src_dir}/{file_path}",
            ]
            print(f"Executing command: {' '.join(command)}")
            subprocess.run(command, check=True)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(
            "Usage: python script.py <src_dir> <path_to_icon_sources> <path_to_integrations.yaml"
        )
        sys.exit(1)

    src_dir = sys.argv[1]
    sdfg_sources_path = sys.argv[2]
    yaml_file_path = sys.argv[3]
    entries = read_yaml(yaml_file_path)
    process_entries(entries, sdfg_sources_path, src_dir)
