from io import StringIO

import requests


def download_and_save_file(url, filename):
    r = requests.get(url, stream=True)

    with open(filename, 'wb') as f:
        for chunk in r.iter_content(chunk_size=1024 * 1024):
            f.write(chunk)


def download_file_stringio(url):
    r = requests.get(url, stream=True)

    file_like = StringIO()
    for chunk in r.iter_content(chunk_size=1024 * 1024, decode_unicode=True):
        file_like.write(chunk)

    file_like.seek(0)
    return file_like


def get_structure_stringio(code):
    return download_file_stringio(f'https://models.rcsb.org/v1/{code}/full')
