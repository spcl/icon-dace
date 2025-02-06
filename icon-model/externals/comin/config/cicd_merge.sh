# installation of "cocogitto" and "pre-commit"

NONINTERACTIVE=1 bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
[ -d /home/linuxbrew/.linuxbrew ] && eval $(/home/linuxbrew/.linuxbrew/bin/brew shellenv)
brew install cocogitto
pip install pre-commit

# perform checks

set -e
cog -v check cocogitto_initial..HEAD
pre-commit run --all-files --config pre-commit-config.merge.yaml
