import re
import requests
from requests.adapters import HTTPAdapter, Retry

re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)


url = 'https://rest.uniprot.org/uniprotkb/search?compressed=true&fields=accession%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clineage%2Clineage_ids%2Cec%2Cannotation_score%2Cgo_p%2Cgo_c%2Cgo%2Cgo_f%2Cprotein_families%2Cxref_eggnog%2Cxref_orthodb%2Cxref_panther%2Cxref_interpro%2Cxref_pfam&format=tsv&query=%28%2A%29&size=500'

print('Downloading annotation data from UNIPROT. This will take few hours.\n')
print('Printing process:\n')
progress = 0

with open('ALL_UNIPROT_ANNOTATION.gz', 'w') as f:
    for batch, total in get_batch(url):
        lines = batch.text.splitlines()
        if not progress:
            print(lines[0], file=f)
        for line in lines[1:]:
            print(line, file=f)
        progress += len(lines[1:])
        print(f'{progress} / {total}')
        
print('Downloading annotation data from UNIPROT finished.')