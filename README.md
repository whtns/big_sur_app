## About
BigSuR is a python dash-based web-application that enables researchers to upload raw scRNA-seq data and perform filtering, analysis, and manual annotation. **It is largely an interactive wrapper for functions provided by [scanpy](https://github.com/theislab/scanpy) and [BigSur](https://github.com/landerlabcode/BigSur)**.

BigSuR is still in a beta state and a few bugs here and there are to be expected. We welcome your bug reports on our issue tracker! 

## Features
BigSuR currently supports the following:

### Data types
* Uploading of 10X mapping directory contents
* Uploading of h5ad formatted scanpy data objects

### Filtering & Feature Selection
* Down-sampling UMI counts
* Filtering by min/max genes per cell
* Filtering by min cells per gene
* Identifying highly-variable genes using **mcFano**, a mean-corrected Fano factor.

### Projections/clustering
* Changing k, the number of neighbors in the neighborhood graph (for projections)
* Changing the clustering resolution used by the [leiden](https://github.com/vtraag/leidenalg) clusting algorithm
* Generating UMAP projections

### Correlation Network Analysis
* Calculation of gene-gene correlation networks using modified Pearson correlation coefficients (mcPCCs).
* Community detection within the correlation graph using the Leiden algorithm.
* Calculation of eigenvector centrality to rank genes within communities.
* Interactive visualization of community centroids and gene-level networks.

### Marker gene analysis
* Detection of marker genes for arbitrary combinations of clusters
* Use of any scanpy-implemented test for marker gene significance

### Annotations
* Viewing clustering and gene expression UMAP projections all at once
* Viewing violin plots of gene expression
* Filtering cells based on any combination of the above plots

### Saving
* Download your anndata object in h5ad format, ready to load into scanpy
for further analysis outside of the bounds of this app

More features are coming! We welcome your suggestions and pull requests.


## Installation

### Linux
There are three main components to the BigSuR software package:
* A redis caching backend server
* A celery task queue for long-running tasks
* The gunicorn server that handles web requests and "runs" the BigSuR server

Redis will need to be installed manually using your distribution's package manager. It does not require any non-default configuration.

We have provided a requirements.txt file that contains a list of all python packages necessary for running BigSuR (including celery and gunicorn). Using pip, install these dependencies: `pip install -r requirements.txt`. 

You will need to manually create the following directory on your computer:
`# mkdir -p /Library/WebServer/Documents/BigSuR/databases`. Then, make sure it is user-writeable using `# chown -R [username_here]:[group_here] /Library/WebServer/Documents/BigSuR`  

Then, run the redis, celery, and gunicorn servers all together using the provided `./BigSuR.sh` script. Point your browser to [the server](http://localhost:8050), and you should be all set to go! 

## Usage tips
* Be patient! Many functions take time (sometimes considerable amounts of it) to process, especially with larger datasets!
* You need not "recalculate everything" when you just want to change the UMAP projection or clustering parameters. Each of those processing sections has a button enabling them to be calculated independently (and thus more rapidly).
* Save your anndata object and load it into scanpy directly for access to any and all of scanpy's extended functionality.

## Development / Debugging
When developing or debugging the web app you may want to run the server in the foreground so interactive debuggers like `ipdb` work as expected.

- Run the app in foreground (recommended for interactive debugging):
```bash
# start redis and celery, then run the web app in foreground
bash ./scripts/start.sh --debug
# or run BigSuR.sh in debug mode
./BigSuR.sh debug
# alternatively run the dev server directly (no gunicorn)
python src/index.py
```

- Run Gunicorn in the foreground for interactive debugging (single sync worker, no timeout):
```bash
python -m gunicorn --pythonpath src -w 1 -k sync -b 0.0.0.0:8050 index:server --timeout 0 --log-level debug
```

- Use a remote debugger (works under Gunicorn without a TTY):
```python
from remote_pdb import RemotePdb
RemotePdb('127.0.0.1', 4444).set_trace()
# then connect from another terminal:
nc 127.0.0.1 4444
```

- Stop the app and cleanup processes safely using the helper script:
```bash
bash scripts/kill_app.sh
# make it executable once if you prefer
chmod +x scripts/kill_app.sh
```

Notes:
- Running the app in the foreground attaches stdin/stdout to your terminal so `ipdb` prompts are interactive.
- By default the provided `start.sh` runs Gunicorn in background; use `--debug` to run the dev server instead.
- If you prefer not to start Celery automatically when debugging, run the dev server directly via `python src/index.py` and start Celery manually if needed.

### License & Origins
This application is based on the **MiCV** app, originally created by [Nigel S. Michki](https://github.com/nigeil) in the [Cai Lab](https://www.cai-lab.org/) at the University of Michigan, and is licensed under the Apache License, Version 2.0. A copy of the license is included in the `LICENSE.txt` file.


This application relies heavily upon the incredible work done by the authors and maintainers of many critical software packages, including:
* [scanpy](https://github.com/theislab/scanpy)
* [BigSur](https://github.com/landerlabcode/BigSur)
* [bbknn](https://github.com/Teichlab/bbknn)
* [louvain/leiden](https://github.com/vtraag/leidenalg)
* [dash](https://plot.ly/dash/)
* [docker](https://www.docker.com/)
* [flask](https://flask.palletsprojects.com)
* [python](https://www.python.org/)
* [numpy](https://numpy.org/)
* [scipy](https://www.scipy.org/)
* [matplotlib](https://matplotlib.org/)
* [seaborn](https://seaborn.pydata.org/)
* And many other dependencies down the line

We thank them for their contributions to open scientific computing and discovery.

## Citations
Please cite the following papers if you use BigSuR in your research:

1. Silkwood K, Dollinger E, Gervin J, Atwood S, Nie Q, Lander AD. Leveraging gene correlations in single cell transcriptomic data. BMC Bioinformatics. 25(1):305. PMCID: PMC11411778
2. Dollinger EP, Silkwood K, Atwood S, Nie Q, Lander AD. Statistically principled feature selection for single cell transcriptomics. BMC Bioinformatics. 26(1):238. PMCID: PMC12490061


