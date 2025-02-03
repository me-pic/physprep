# Contribution

## Cloning the project from Github

1. Fork the repository and clone the project:

```bash
git clone git@github.com:<your_username>/physprep.git
```

2. Create a virtual environment:

```bash
python3 -m venv <venv_name>
source <venv_name>/bin/activate
```

3. Install the dependency:

```bash
cd physprep
pip install -e .[dev]
```

4. Install pre-commit hooks
```bash
pre-commit install
```