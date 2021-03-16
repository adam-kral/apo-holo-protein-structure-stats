import concurrent.futures
import logging
import urllib.request

import pandas as pd

from apo_holo_structure_stats.paper_repl.main import download_structure

if __name__ == '__main__':
    logging.root.setLevel(logging.INFO)

    df = pd.read_csv('holo.dat', delimiter=r'  ', comment='#', header=None,
                     names=('domain_count', 'holos_count', 'holos'))

    pdb_codes = [paper_code[:4] for index, row in df.iterrows() for paper_code in row.holos.split()]

    opener = urllib.request.build_opener()
    opener.addheaders = [('User-agent', 'Mozilla/5.0 (Macintosh; Intel Mac OS X 11.2; rv:86.0) Gecko/20100101 Firefox/86.0')]

    with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:  # >= 1 workers 421 too many connections, dneska už nic nestáhnu..
       executor.map(download_structure, pdb_codes)
