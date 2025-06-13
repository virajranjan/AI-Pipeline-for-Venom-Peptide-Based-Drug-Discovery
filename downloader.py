import requests
import os


BASE_URL = "https://rest.uniprot.org/uniprotkb/stream"
QUERY = 'keyword:"Toxin" AND taxonomy_id:6493'
OUTPUT_FILE = "/home/zorro/project_peptide/fasta/conus_toxins_full_uniprot.fasta"

def download_all_fasta():
    headers = {"User-Agent": "ConoVenomDownloader/1.0"}
    params = {
        "query": QUERY,
        "format": "fasta"
    }

    print("[INFO] Requesting all conotoxin sequences from UniProt...")
    response = requests.get(BASE_URL, params=params, headers=headers)

    if response.status_code != 200:
        print(f"[ERROR] Failed with status code {response.status_code}")
        print("[DETAIL]", response.text)
        return

    with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
        f.write(response.text)

    num_sequences = response.text.count(">")
    print(f"[DONE] Downloaded {num_sequences} FASTA sequences.")
    print(f"[SAVED] Output written to: {OUTPUT_FILE}")


if __name__ == "__main__":
    download_all_fasta()
