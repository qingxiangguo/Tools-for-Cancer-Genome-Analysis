# Poetry: A Practical Guide

Poetry is a Python dependency management and packaging tool, designed to be user-friendly and intuitive, while ensuring consistency and reproducibility in projects.

## Installation

Before diving into the commands, ensure you've installed Poetry:

```bash
curl -sSL https://install.python-poetry.org | python3 -
```

Or using pip (though the above method is recommended):

```bash
pip install poetry
```

## Key Concepts and Commands

### 1. Starting a New Project

```bash
poetry new project-name
```

This command creates a new directory `project-name` with a basic project structure.

### 2. Adding Dependencies

```bash
poetry add package-name
```

This command adds a package to your project. Poetry will automatically update the `pyproject.toml` and `poetry.lock` files.

For development-only dependencies (like testing libraries):

```bash
poetry add package-name --dev
```

### 3. Removing Dependencies

```bash
poetry remove package-name
```

This will remove the package and update the necessary files.

### 4. Installing Dependencies

If you have an existing `pyproject.toml` and `poetry.lock` file (e.g., after cloning a repo):

```bash
poetry install
```

This command reads the `poetry.lock` file to ensure consistent installations across environments.

### 5. Building and Publishing

To build your project:

```bash
poetry build
```

This produces distribution packages (wheel and sdist) that can be published to PyPI.

To publish to PyPI:

```bash
poetry publish --build
```

Ensure you're authenticated to PyPI before publishing.

### 6. Running Scripts

You can define custom scripts in `pyproject.toml`:

```toml
[tool.poetry.scripts]
script-name = "module:function"
```

Run the script using, this will run in poetry Virtual Environment

```bash
poetry run script-name
```

### 7. Updating Dependencies

To update a specific package and simultaneously refresh both the `poetry.lock` and `pyproject.toml` files:

```bash
poetry update package-name
```

To update all packages and refresh the `poetry.lock`:

```bash
poetry update
```

#### Handling the `poetry.lock` File

If you've manually updated the `pyproject.toml` but not the `poetry.lock`:

1. Delete `poetry.lock` and then run:

```bash
poetry install
```

This will regenerate the `poetry.lock` file.

2. Alternatively, you can refresh all packages and the `poetry.lock` with:

```bash
poetry update
```

3. If you want to regenerate the lock file without updating packages:

```bash
poetry lock --no-update
```

### 8. Show Dependencies

To list all dependencies and their details:

```bash
poetry show
```

### 9. Virtual Environments # So you don't need to run: poetry run python XXX every-time.

Poetry creates a virtual environment for your projects by default. You can access the environment using:

```bash
poetry shell
```

### 10. Check

To check if `pyproject.toml` is valid:

```bash
poetry check
```

## Conclusion

Poetry is a robust tool that makes dependency management and packaging seamless. With its deterministic builds, it ensures that projects are consistent across all environments. Get started, and enjoy a better Python packaging experience!
