import requests, json, time

all_records = []
limit = 1000
offset = 0
select = "structure_id,linear_sequence,linear_sequence_length,curated_source_antigens,qualitative_measures,antibody_isotypes,mhc_allele_names,host_organism_names,reference_ids"

while True:
    t0 = time.time()
    r = requests.get("https://query-api.iedb.org/epitope_search",
                      params={"source_organism_iri_search": "cs.{NCBITaxon:12637}",
                              "structure_type": "eq.Linear peptide",
                              "order": "structure_id",
                              "select": select,
                              "limit": limit, "offset": offset},
                      timeout=120)
    r.raise_for_status()
    batch = r.json()
    all_records.extend(batch)
    print(f"offset={offset} got={len(batch)} total_so_far={len(all_records)} took={time.time()-t0:.1f}s", flush=True)
    if len(batch) < limit:
        break
    offset += limit
    time.sleep(0.3)

with open("references/iedb_dengue_linear_epitopes.json", "w", encoding="utf-8") as f:
    json.dump(all_records, f)

print("TOTAL:", len(all_records))
