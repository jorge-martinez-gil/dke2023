# Contributing to Neurofuzzy

Thank you for contributing to this repository.

## Development setup

1. Create and activate the environment:
   ```bash
   conda env create -f environment.yml
   conda activate neurofuzzy
   ```
2. Install in editable mode with development tools:
   ```bash
   pip install -e ".[dev]"
   ```

## Running tests

Run the full test suite with:

```bash
pytest tests/ -v
```

## Code style

- Follow PEP 8.
- Use type hints in public APIs.
- Add NumPy-style docstrings for modules, classes, and functions.

## Adding a new dataset

1. Place the dataset file in `datasets/`.
2. Use the same format: first column ground truth, remaining columns features.
3. Add the dataset key/path in `experiments/run_experiment.py`.
4. Add tests validating shape and load behavior.

## Adding a new defuzzification method

1. Extend `neurofuzzy/fuzzy_controller.py` with the new method implementation.
2. Add an argument or selector in `FuzzyController` and `NeurofuzzyModel`.
3. Add tests that validate numerical outputs and output range `[0, 1]`.
4. Document the method and assumptions in `README.md`.

## Pull request checklist

- [ ] Tests pass locally (`pytest tests/ -v`)
- [ ] New code includes type hints and NumPy-style docstrings
- [ ] README/Docs updated if behavior changed
- [ ] No generated artifacts committed (`results/`, `*.npy`, caches)
