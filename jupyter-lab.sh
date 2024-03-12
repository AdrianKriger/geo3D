#!/usr/bin/env bash
echo "🪛 Running Jupyter using Venv:"
echo "Switch the kernel to 'virtualenv'"
echo "after opening in its own window"
echo "--------------------------------"
python -m ipykernel install --user --name=virtualenv
jupyter kernelspec list
jupyter-lab
