#!/usr/bin/env python3
"""
build_mouse_gse_hdf5.py

Reconstructs the GEO MINiML metadata dump that filter_gse_minimal_metadata.py
consumes, scoped to mouse instead of human.

Produces an HDF5 with a single dataset 'GSE_XML_strings': an indexed array where
each entry is the full family MINiML XML (as a UTF-8 string) for one GSE. This is
the exact layout the downstream filtering script reads via
    in_file['GSE_XML_strings']
so no changes to that script's HDF5 access are needed.

Two stages:
  1. eSearch the 'gds' Entrez database for every mouse GSE under
     "expression profiling by high throughput sequencing", page through all UIDs,
     eSummary them to GSE accessions.
  2. For each GSE, download its family MINiML from the GEO FTP tree, extract the
     XML, and write the full XML string into the HDF5.

Run this on Oscar (or anywhere with network + h5py). It is an overnight job at
mouse scale (tens of thousands of series). Checkpointing is built in so you can
resume: stage 1 writes an accession list; stage 2 skips GSEs already written.

Usage:
    # Set NCBI_API_KEY in your environment first (raises eSearch rate limit 3->10/s)
    export NCBI_API_KEY=your_key_here
    export NCBI_EMAIL=you@brown.edu
    python build_mouse_gse_hdf5.py

Dependencies: requests, h5py  (both already in order_env, or pip install)
"""

import os
import io
import gzip
import time
import tarfile
import urllib.request
import urllib.error
from os.path import join, exists
from xml.etree import ElementTree as ET

import requests
import h5py

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
DATA_DIR = '/users/dhan30/splicing_order/data'

# Stage 1 output: the list of mouse GSE accessions matching the query
ACCESSION_LIST = join(DATA_DIR, 'mouse_expression_profiling_by_HT_seq_gse_accessions.txt')

# Stage 2 output: the HDF5 dump the filtering script reads
HDF5_OUT = join(DATA_DIR, 'mouse_expression_profiling_by_HT_seq_gds_metadata.hdf5')

# The GEO DataSets query. This is the mouse analogue of the human dump's query.
# [ETYP] restricts DataSet type; gse[Entry Type] restricts results to Series
# (not GDS/GPL/GSM) so each UID maps to one GSE.
ENTREZ_QUERY = (
    '"Mus musculus"[Organism] AND '
    '"expression profiling by high throughput sequencing"[DataSet Type] AND '
    'gse[Entry Type]'
)

EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
API_KEY = os.environ.get('NCBI_API_KEY', '')
EMAIL = os.environ.get('NCBI_EMAIL', '')
TOOL = 'splicing_order_mouse_build'

# Politeness: with an API key NCBI allows 10 req/s; without, 3 req/s.
REQ_DELAY = 0.11 if API_KEY else 0.34
RETMAX = 500  # UIDs per eSearch page

XMLNS = '{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}'


def _eutils_params(**extra):
    p = {'tool': TOOL}
    if API_KEY:
        p['api_key'] = API_KEY
    if EMAIL:
        p['email'] = EMAIL
    p.update(extra)
    return p


# ---------------------------------------------------------------------------
# Stage 1: query GEO DataSets -> list of GSE accessions
# ---------------------------------------------------------------------------
def stage1_get_gse_accessions():
    if exists(ACCESSION_LIST):
        with open(ACCESSION_LIST) as f:
            accs = [line.strip() for line in f if line.strip()]
        print(f'[stage1] Reusing existing accession list: {len(accs)} GSEs')
        return accs

    print('[stage1] eSearch gds for mouse HT-seq series...')
    # Use history server so we can page through all results.
    r = requests.get(
        f'{EUTILS}/esearch.fcgi',
        params=_eutils_params(db='gds', term=ENTREZ_QUERY,
                              usehistory='y', retmax=0),
        timeout=60,
    )
    r.raise_for_status()
    root = ET.fromstring(r.text)
    count = int(root.findtext('Count'))
    webenv = root.findtext('WebEnv')
    query_key = root.findtext('QueryKey')
    print(f'[stage1] {count} total results to page through')

    # eSummary straight off the history server: pass query_key + WebEnv + retstart
    # instead of a UID list. This avoids building a multi-kilobyte id=... URL
    # (the source of the 414 error) entirely. We POST as belt-and-suspenders,
    # since NCBI recommends POST for any request carrying many identifiers.
    accessions = []
    for retstart in range(0, count, RETMAX):
        time.sleep(REQ_DELAY)
        r = requests.post(
            f'{EUTILS}/esummary.fcgi',
            data=_eutils_params(db='gds', query_key=query_key, WebEnv=webenv,
                                retstart=retstart, retmax=RETMAX),
            timeout=120,
        )
        r.raise_for_status()
        summ = ET.fromstring(r.text)
        for docsum in summ.findall('.//DocSum'):
            acc = None
            entry_type = None
            for item in docsum.findall('Item'):
                if item.attrib.get('Name') == 'Accession':
                    acc = item.text
                elif item.attrib.get('Name') == 'entryType':
                    entry_type = item.text
            # Keep only GSE-level records
            if acc and acc.startswith('GSE') and (entry_type in (None, 'GSE')):
                accessions.append(acc)
        print(f'[stage1]   summarized {min(retstart + RETMAX, count)}/{count}')

    # De-duplicate, preserve order
    seen = set()
    accessions = [a for a in accessions if not (a in seen or seen.add(a))]

    os.makedirs(DATA_DIR, exist_ok=True)
    with open(ACCESSION_LIST, 'w') as f:
        for a in accessions:
            f.write(a + '\n')
    print(f'[stage1] Wrote {len(accessions)} unique GSE accessions -> {ACCESSION_LIST}')
    return accessions


# ---------------------------------------------------------------------------
# Stage 2: per-GSE MINiML fetch -> HDF5
# ---------------------------------------------------------------------------
def _miniml_ftp_url(gse):
    # GSE12345 -> GSE12nnn ; GSE1 -> GSEnnn
    num = gse[3:]
    if len(num) <= 3:
        stub = 'GSEnnn'
    else:
        stub = 'GSE' + num[:-3] + 'nnn'
    return (f'https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/'
            f'miniml/{gse}_family.xml.tgz')


def _fetch_miniml_xml(gse, retries=3):
    """Download <GSE>_family.xml.tgz, return the family XML as a UTF-8 string."""
    url = _miniml_ftp_url(gse)
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=120) as resp:
                raw = resp.read()
            tar = tarfile.open(fileobj=io.BytesIO(raw), mode='r:gz')
            # The family XML is the member ending in _family.xml
            member = next((m for m in tar.getmembers()
                           if m.name.endswith('_family.xml')), None)
            if member is None:
                return ''
            xml_bytes = tar.extractfile(member).read()
            return xml_bytes.decode('utf-8', errors='replace')
        except (urllib.error.URLError, tarfile.TarError, OSError) as e:
            if attempt == retries - 1:
                print(f'[stage2]   WARN {gse}: {e} (giving up, writing empty)')
                return ''
            time.sleep(2 ** attempt)
    return ''


def stage2_build_hdf5(accessions):
    n = len(accessions)

    # Resume support: if the HDF5 exists, find which slots are already filled.
    start_idx = 0
    if exists(HDF5_OUT):
        with h5py.File(HDF5_OUT, 'r') as f:
            if 'GSE_XML_strings' in f:
                filled = f['GSE_XML_strings']
                # First empty slot = resume point
                for i in range(len(filled)):
                    if len(filled[i]) == 0:
                        start_idx = i
                        break
                else:
                    start_idx = len(filled)
        print(f'[stage2] Resuming from index {start_idx}/{n}')

    str_dt = h5py.string_dtype(encoding='utf-8')
    mode = 'a' if exists(HDF5_OUT) else 'w'
    with h5py.File(HDF5_OUT, mode) as f:
        if 'GSE_XML_strings' not in f:
            ds = f.create_dataset('GSE_XML_strings', shape=(n,),
                                  dtype=str_dt, compression='gzip')
        else:
            ds = f['GSE_XML_strings']

        for i in range(start_idx, n):
            gse = accessions[i]
            time.sleep(REQ_DELAY)
            ds[i] = _fetch_miniml_xml(gse)
            if (i + 1) % 50 == 0:
                f.flush()
                print(f'[stage2]   {i + 1}/{n} series written')
        f.flush()
    print(f'[stage2] Done. HDF5 -> {HDF5_OUT}')


if __name__ == '__main__':
    accs = stage1_get_gse_accessions()
    stage2_build_hdf5(accs)