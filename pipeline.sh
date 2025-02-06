SCRATCHPAD_DIR="icon-scratchpad"
SOURCE_DIR="icon-model"

# Remove the icon-scratchpad directory only if it exists
if [ -d "$SCRATCHPAD_DIR" ]; then
    rm -r "$SCRATCHPAD_DIR"
fi

# Check if the directory exists and is not empty
if [ -d "$SCRATCHPAD_DIR" ] && [ "$(ls -A $SCRATCHPAD_DIR)" ]; then
    echo "Error: Directory '$SCRATCHPAD_DIR' exists and is not empty."
    exit 1
fi

mkdir -p "$SCRATCHPAD_DIR"

echo "Directory '$SCRATCHPAD_DIR' created or already exists and is empty."

ABS_SOURCE_DIR=$(realpath "$SOURCE_DIR")
ABS_SCRATCHPAD_DIR=$(realpath "$SCRATCHPAD_DIR")

cp -sR "$ABS_SOURCE_DIR" "$ABS_SCRATCHPAD_DIR"

# Run folder should not be symlink
rm -R "$ABS_SCRATCHPAD_DIR/icon-model/run"
cp -R "$ABS_SOURCE_DIR/run" "$ABS_SCRATCHPAD_DIR/icon-model"


echo "Contents of 'icon-model' copied to '$SCRATCHPAD_DIR' as symlinks."

python generate_sdfgs.py "$ABS_SCRATCHPAD_DIR/$SOURCE_DIR" "$ABS_SCRATCHPAD_DIR/$SOURCE_DIR/sdfgs" "$ABS_SCRATCHPAD_DIR/$SOURCE_DIR/sdfgs/integrations.yaml"

#python generate.py get_albedos ./integrations.yaml .
#python optimize_get_albedos.py ./get_albedos_unsimplified.sdfgz ./get_albedos_optimized.sdfgz