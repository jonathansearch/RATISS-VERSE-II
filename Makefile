.PHONY: build verify test clean

build:
	python scripts/build.py

verify:
	python scripts/verify.py

test:
	python -m pytest tests/ -v

clean:
	rm -f results/*.json results/*.png figures/*.png results/MANIFEST.sha256
	rm -rf .pytest_cache __pycache__ src/v2/__pycache__ tests/__pycache__
