# Notebooks

| Notebook | Audience | Contents |
|---|---|---|
| [`01_getting_started.ipynb`](01_getting_started.ipynb) | Beginner | Load data, fit the model, evaluate, compare to a baseline, visualize, and observe the train/test gap. Includes exercises. |

Run with Jupyter from this directory:

```bash
pip install -e ".[viz]" jupyter
jupyter lab notebooks/01_getting_started.ipynb
```

The notebook uses relative paths (`../datasets/...`) and a small optimizer budget
so it completes in a couple of minutes.
