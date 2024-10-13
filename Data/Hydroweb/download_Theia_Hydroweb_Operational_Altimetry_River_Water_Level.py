help_message = """
Download products from your Hydroweb.next project (https://hydroweb.next.theia-land.fr) using EODAG (https://github.com/CS-SI/eodag)
This script is an example tuned for your last Hydroweb.next project but feel free to adapt it for future requests.
Follow these steps:
1. If not already done, install EODAG and packaging latest version using `pip install -U eodag packaging` or `conda update eodag packaging`
2a. Generate an API-Key from Hydroweb.next portal in your user settings
2b. Carefully store your API-Key
- either in your eodag configuration file (usually ~/.config/eodag/eodag.yml, automatically generated the first time you use eodag) in auth/credentials/apikey="PLEASE_CHANGE_ME"
- or in an environment variable `export EODAG__HYDROWEB_NEXT__AUTH__CREDENTIALS__APIKEY="PLEASE_CHANGE_ME"`
3. You can change download directory by modifying the variable path_out. By default, current path is used.
4. You are all set, run this script `python download_Theia_Hydroweb_Operational_Altimetry_River_Water_Level.py`

For more information, please refer to EODAG Documentation https://eodag.readthedocs.io
"""
try:
    import os
    import logging
    from eodag import EODataAccessGateway, SearchResult
    from eodag import __version__ as eodag_version
    from eodag import setup_logging
    from packaging import version

    assert version.parse(eodag_version) >= version.parse("2.11.0"), help_message
except ImportError:
    print(help_message)
    exit(1)

setup_logging(1)  # 0: nothing, 1: only progress bars, 2: INFO, 3: DEBUG

# Add custom logger
custom_logger = logging.getLogger("Hydroweb.next")
custom_logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)-15s %(name)-32s [%(levelname)-8s] %(message)s')
handler = logging.StreamHandler()  # Use a different handler if needed
handler.setFormatter(formatter)
custom_logger.addHandler(handler)

# Set timeout to 30s
os.environ["EODAG__HYDROWEB_NEXT__SEARCH__TIMEOUT"] = "30"

dag = EODataAccessGateway()

# --------------------------------------------------
# Declare the path where the file will be downloaded
path_out = f"{os.getcwd()}/results"

if not os.path.isdir(path_out):
    raise RuntimeError(f"Path {path_out} does not exist")

# Default search criteria when iterating over collection pages
default_search_criteria = {
    "items_per_page": 2000,
}

# ---------------------------------------------------------------------------------------------------------------------
# This command will perform a search using provided query arguments.
# - specify a collection using the `productType` key
# - add time restrictions using the `start` and `end` keys (e.g. "start": "2020-05-01" , "end": "2020-05-10T00:00:00Z",
#   UTC ISO8601 format)
# - add spatial restrictions using the "geom" key (e.g. "geom": "POLYGON ((1 43, 2 43, 2 44, 1 44, 1 43))" WKT string,
#   a bounding-box list [lonmin, latmin, lonmax, latmax] can also be passed )
# - more query arguments can be used, see
#   https://eodag.readthedocs.io/en/stable/notebooks/api_user_guide/4_search.html?#Search-parameters
query_args = {
    "productType": "HYDROWEB_RIVERS_OPE",
    "geom": [72, 20, 100, 33]
}
# Niger: [-15, 2, 15, 20]
# Ganges-Brahmaputra: [72, 20, 100, 33]

query_args.update(default_search_criteria)

# Iterate over all pages to find all products
search_results = ([])
custom_logger.info(f"Searching for products matching {query_args}...")

for i, page_results in enumerate(dag.search_iter_page(**query_args)):
    custom_logger.info(f"{len(page_results)} product(s) found on page {i+1}")
    search_results.extend(page_results)

custom_logger.info(f"Total products found : {len(search_results)}")

# Download all found products
custom_logger.info(f"Downloading {len(search_results)} products...")

# This command actually downloads the matching products
downloaded_paths = dag.download_all(search_results, outputs_prefix=path_out)
if downloaded_paths:
    distinct_values = list(set(downloaded_paths))
    # Check is distinct values length is equal to all results length
    if len(distinct_values) != len(search_results):
        custom_logger.warning(f"Distinct values length is not equal to all results length. {len(distinct_values)} != {len(search_results)}")
        custom_logger.warning(f"Some files have not been downloaded")
    else:
        custom_logger.info(f"All {len(search_results)} files have been successfully downloaded.")
else:
    print(f"No files downloaded! Verify API-KEY and/or product search configuration.")

