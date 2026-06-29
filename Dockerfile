# Reproducible environment for the neurofuzzy benchmark.
# Build:  docker build -t neurofuzzy .
# Run all:  docker run --rm -v "$PWD/benchmarks:/app/benchmarks" neurofuzzy
FROM python:3.10-slim

ENV PYTHONUNBUFFERED=1 \
    MPLCONFIGDIR=/tmp/mplconfig \
    PIP_NO_CACHE_DIR=1

WORKDIR /app

COPY pyproject.toml requirements-lock.txt README.md ./
COPY neurofuzzy ./neurofuzzy
COPY datasets ./datasets
COPY experiments ./experiments
COPY tests ./tests

RUN pip install -r requirements-lock.txt && pip install -e ".[viz]"

# Default: run the full benchmark and regenerate figures into mounted volume.
CMD ["sh", "-c", "python -m neurofuzzy.benchmark --all --seeds 10 --out benchmarks/results && python -m neurofuzzy.visualize --results benchmarks/results --out benchmarks/figures"]
