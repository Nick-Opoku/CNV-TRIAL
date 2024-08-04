import csv

def parse_refFlat(file_path):
    refFlat_records = []
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            gene_name, transcript_id, chromosome, strand, tx_start, tx_end = row[0], row[1], row[2], row[3], int(row[4]), int(row[5])
            refFlat_records.append({
                "gene_name": gene_name,
                "transcript_id": transcript_id,
                "chromosome": chromosome,
                "strand": strand,
                "tx_start": tx_start,
                "tx_end": tx_end
            })
    return refFlat_records

def parse_cnv(file_path):
    cnv_records = []
    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip header
        for row in reader:
            chromosome, st_bp, end_bp, length_kb, cnv_type = int(row[0]), int(row[1]), int(row[2]), float(row[3]), row[4]
            cnv_records.append({
                "chromosome": chromosome,
                "st_bp": st_bp,
                "end_bp": end_bp,
                "length_kb": length_kb,
                "cnv_type": cnv_type
            })
    return cnv_records

def find_overlapping_genes(cnv_records, refFlat_records):
    results = []
    for cnv in cnv_records:
        cnv_chromosome = f"chromosome {cnv['chromosome']}"
        cnv_start, cnv_end = cnv['st_bp'], cnv['end_bp']
        overlapping_genes = []
        for gene in refFlat_records:
            gene_chromosome = gene['chromosome']
            gene_start, gene_end = gene['tx_start'], gene['tx_end']
            if cnv_chromosome == gene_chromosome and not (gene_end < cnv_start or gene_start > cnv_end):
                overlapping_genes.append(gene['gene_name'])
        results.append({
            "chromosome": cnv['chromosome'],
            "st_bp": cnv['st_bp'],
            "end_bp": cnv['end_bp'],
            "length_kb": cnv['length_kb'],
            "cnv_type": cnv['cnv_type'],
            "overlapping_genes": overlapping_genes
        })
    return results

def write_results_to_csv(results, output_file):
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['chromosome', 'st_bp', 'end_bp', 'length_kb', 'cnv_type', 'overlapping_genes']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        for result in results:
            writer.writerow({
                'chromosome': result['chromosome'],
                'st_bp': result['st_bp'],
                'end_bp': result['end_bp'],
                'length_kb': result['length_kb'],
                'cnv_type': result['cnv_type'],
                'overlapping_genes': ", ".join(result['overlapping_genes'])
            })

def main():
    refFlat_file = 'sorted_refFlat.txt'
    cnv_file = 'resultsFIN.csv'
    output_file = 'merged.csv'

    refFlat_records = parse_refFlat(refFlat_file)
    cnv_records = parse_cnv(cnv_file)
    overlapping_genes_results = find_overlapping_genes(cnv_records, refFlat_records)

    # Write results to CSV
    write_results_to_csv(overlapping_genes_results, output_file)

if __name__ == "__main__":
    main()

