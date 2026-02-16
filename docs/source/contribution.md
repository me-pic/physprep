# Contribution

## Cloning the project from Github

1. Fork the repository and clone the project:

```bash
git clone git@github.com:<your_username>/physprep.git
```

## Installing the dependencies

### Using poetry

1. Create a virtual environment and install the dependencies:

```bash
poetry install --with dev
```
2. Activate the environment

```bash
poetry shell
```

If the command above does not work, use:

```bash
source $(poetry env info --path)/bin/activate
```

### Using venv

1. Create a virtual environment:

```bash
python3 -m venv <venv_name>
source <venv_name>/bin/activate
```

2. Install the dependency:

```bash
cd physprep
pip install -e .[dev]
```

3. Install pre-commit hooks
```bash
pre-commit install
```

## Building the documentation locally

In your terminal:

```bash
cd docs
make html
```

The documentation will be generated in `docs/build/html`.
