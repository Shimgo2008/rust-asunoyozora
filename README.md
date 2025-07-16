# source
https://github.com/Shimgo2008/rosic-asunoyozora/

# build command
build binary
```sh
cargo build --release
```

build pymodule
```sh
python -m venv .venv
source .venv/bin/activate
pip install -r requirement.txt
maturin build
```