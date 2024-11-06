import requests

def get_snp_rs(chromosome, start, end):
    base_url_search = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
    base_url_summary = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
    
    # NCBI parameters for SNP search
    search_params = {
        "db": "snp",
        "term": f"{chromosome}[CHR] AND {start}:{end}[CHRPOS]",
        "retmode": "json",
    }
    
    # Send request to NCBI to search for SNPs in the given region
    search_response = requests.get(base_url_search, params=search_params)
    search_data = search_response.json()
    
    # Extract the list of SNP IDs (rs numbers)
    snp_ids = search_data['esearchresult']['idlist']
    
    # Fetch detailed information for each SNP ID
    if not snp_ids:
        return None
    
    summary_params = {
        "db": "snp",
        "id": ",".join(snp_ids),
        "retmode": "json"
    }
    summary_response = requests.get(base_url_summary, params=summary_params)
    summary_data = summary_response.json()
    
    # Look for the canonical SNP (non-merged)
    for snp_id in snp_ids:
        snp_info = summary_data['result'][snp_id]
        if 'mergedrs' not in snp_info:
            return snp_info['snp_id']  # Return the canonical SNP rs number
    
    return None  # No canonical SNP found

# Example usage
chromosome = 1
start = 841085
end = 841085

snp_rs_id = get_snp_rs(chromosome, start, end)
print("Main RS Number:", snp_rs_id)






#%%
import requests

def get_snp_rs(chromosome, start, end):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
    
    # NCBI parameters for SNP search
    params = {
        "db": "snp",
        "term": f"{chromosome}[CHR] AND {start}:{end}[CHRPOS]",
        "retmode": "json",
    }
    
    # Send request to NCBI to search for SNPs in the given region
    response = requests.get(base_url, params=params)
    data = response.json()
    
    # Extract the list of SNP IDs (rs numbers)
    snp_ids = data['esearchresult']['idlist']
    
    # Return the SNP rs numbers
    return snp_ids

# Example usage
chromosome = 1
start = 841085
end = 841085

snp_rs_ids = get_snp_rs(chromosome, start, end)
print("RS Numbers:", snp_rs_ids)




#%%
import requests

def get_snp_rs(chromosome, start, end):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
    
    # NCBI parameters for SNP search
    params = {
        "db": "snp",
        "term": f"{chromosome}[CHR] AND {start}:{end}[CHRPOS]",
        "retmode": "json",
    }
    
    # Send request to NCBI to search for SNPs in the given region
    response = requests.get(base_url, params=params)
    data = response.json()
    
    # Extract the list of SNP IDs (rs numbers)
    snp_ids = data['esearchresult']['idlist']
    
    return snp_ids

def get_snp_details(snp_ids):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
    
    # Prepare the parameters for the SNP details request
    params = {
        "db": "snp",
        "id": ",".join(snp_ids),  # Join SNP IDs as comma-separated string
        "retmode": "json"
    }
    
    # Send request to NCBI to get SNP details
    response = requests.get(base_url, params=params)
    data = response.json()
    
    # Extract details for each SNP
    snp_details = []
    for snp_id in snp_ids:
        if snp_id in data['result']:
            snp_info = data['result'][snp_id]
            snp_details.append({
                "rs_number": snp_info.get("uid"),
                "variant_type": snp_info.get("snp_class"),
                "alleles": snp_info.get("global_molecule"),
                "gene": snp_info.get("gene", {}).get("name", "N/A")
            })
    
    return snp_details

# Example usage
chromosome = 1
start = 70352
end = 70352

# Step 1: Get the SNP rs numbers in the given region
snp_rs_ids = get_snp_rs(chromosome, start, end)

# Step 2: Get detailed information about each SNP
snp_details = get_snp_details(snp_rs_ids)

# Print the SNP details
for snp in snp_details:
    print(f"RS Number: {snp['rs_number']}")
    print(f"Variant Type: {snp['variant_type']}")
    print(f"Alleles: {snp['alleles']}")
    print(f"Gene: {snp['gene']}")
    print("-" * 30)



#%%
import requests

def get_snp_rs(chromosome, start, end):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
    
    # NCBI parameters for SNP search
    params = {
        "db": "snp",
        "term": f"{chromosome}[CHR] AND {start}:{end}[CHRPOS]",
        "retmode": "json",
    }
    
    # Send request to NCBI to search for SNPs in the given region
    response = requests.get(base_url, params=params)
    data = response.json()
    
    # Extract the list of SNP IDs (rs numbers)
    snp_ids = data['esearchresult']['idlist']
    
    return snp_ids

def get_snp_details(snp_ids):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
    
    # Prepare the parameters for the SNP details request
    params = {
        "db": "snp",
        "id": ",".join(snp_ids),  # Join SNP IDs as comma-separated string
        "retmode": "json"
    }
    
    # Send request to NCBI to get SNP details
    response = requests.get(base_url, params=params)
    data = response.json()
    
    # Extract details for each SNP
    snp_details = []
    for snp_id in snp_ids:
        if snp_id in data['result']:
            snp_info = data['result'][snp_id]
            
            # Attempt to extract alleles and gene information from the relevant fields
            alleles = snp_info.get("global_molecule")
            gene = snp_info.get("gene", {}).get("name") if snp_info.get("gene") else snp_info.get("genes", [{}])[0].get("name", "N/A")
            snp_details.append({
                "rs_number": snp_info.get("uid"),
                "variant_type": snp_info.get("snp_class"),
                "alleles": alleles if alleles else "Not available",
                "gene": gene if gene else "N/A"
            })
    
    return snp_details

# Example usage
chromosome = 1
start = 70352
end = 70352

# Step 1: Get the SNP rs numbers in the given region
snp_rs_ids = get_snp_rs(chromosome, start, end)

# Step 2: Get detailed information about each SNP
snp_details = get_snp_details(snp_rs_ids)

# Print the SNP details
for snp in snp_details:
    print(f"RS Number: {snp['rs_number']}")
    print(f"Variant Type: {snp['variant_type']}")
    print(f"Alleles: {snp['alleles']}")
    print(f"Gene: {snp['gene']}")
    print("-" * 30)






















#%%
import pandas as pd
import requests
import time
from requests.exceptions import ConnectionError, Timeout

# Function to get rs number based on chr, pos, ref, alt
def get_rs_number(chr, pos, ref, alt, assembly='GRCh38', retries=5, delay=1):
    # Check if ref or alt alleles are more than one letter
    #if len(ref) > 1 or len(alt) > 1:
    #    return None
    
    url = f"http://rest.ensembl.org/overlap/region/human/{chr}:{pos}-{pos}?feature=variation;content-type=application/json;assembly={assembly}"
    
    for attempt in range(retries):
        try:
            response = requests.get(url, headers={"Content-Type": "application/json"}, timeout=10)
            if response.ok:
                data = response.json()
                
                # Extract rsID if available and check alleles
                if data:
                    for variant in data:
                        if 'id' in variant and variant['id'].startswith('rs'):
                            alleles = variant.get('alleles', [])
                            alleles_list = [allele.split('/')[0] for allele in alleles]
                            
                            # Check if both ref and alt are in the found alleles
                            if ref in alleles_list and alt in alleles_list:
                                return variant['id']
            else:
                response.raise_for_status()
        except (ConnectionError, Timeout) as e:
            print(f"Attempt {attempt+1} failed with error: {e}. Retrying in {delay} seconds...")
            time.sleep(delay)
    return None

# Example usage
chr = '1'
pos = 113687822
ref = 'AC'
alt = 'A'
rs_id = get_rs_number(chr, pos, ref, alt)
print(f"RS Number for {chr}:{pos} with alleles {ref}/{alt} is {rs_id}")

# Initialize an empty DataFrame to store the results for all chromosomes
input_data = pd.DataFrame()

# Loop over chromosomes 1 to 22
for chr_num in range(1, 23):
    # Define the file path for the GWAS data for each chromosome
    data_file = f"/Users/fvashegh/Library/CloudStorage/OneDrive-JNJ/Pi/GWAS/UKB.bt.htp.IMP.chr{chr_num}.all.additive.SNV.28_X714.regenie"
    
    # Read the GWAS data
    data = pd.read_csv(data_file, sep="\t")
    
    # Filter by p-value
    filtered_data = data[data['Pval'] <= 5e-5]
    
    # Select only the specified columns
    filtered_data = filtered_data[['Chr', 'Pos', 'Ref', 'Alt', 'Pval']]
    
    # Remove rows where 'Ref' or 'Alt' have values with length greater than 1 (e.g. 'TCTC')
    filtered_data = filtered_data[filtered_data['Ref'].str.len() == 1]
    filtered_data = filtered_data[filtered_data['Alt'].str.len() == 1]
    
    # Append the results to the final DataFrame
    input_data = pd.concat([input_data, filtered_data], ignore_index=True)

# Apply the get_rs_number function to each row and add the rs number as a new column
input_data['snp'] = input_data.apply(
    lambda row: get_rs_number(row['Chr'], row['Pos'], row['Ref'], row['Alt']), axis=1
)

# Rename columns to lowercase
input_data.rename(columns={
    'Chr': 'chr',
    'Pos': 'pos',
    'Ref': 'ref',
    'Alt': 'alt',
    'Pval': 'pvalue'
}, inplace=True)

# Reorder columns to bring 'snp' as the first column
input_data = input_data[['snp', 'chr', 'pos', 'ref', 'alt', 'pvalue']]

# Display the first few rows of the final data with the new 'snp' column as the first
print(input_data.head())

# Save the input_data DataFrame as a CSV file
csv_file_path = '/Users/fvashegh/Library/CloudStorage/OneDrive-JNJ/Pi/GWAS/RA_GWAS.csv'
input_data.to_csv(csv_file_path, index=False)



#%%
import requests
import urllib3
import time
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
import pyreadr

# Disable SSL warnings
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

def get_snp_rs_grch37(chromosome, position, retries=10, delay=1):
    base_url = f"https://grch37.rest.ensembl.org/overlap/region/human/{chromosome}:{position}-{position}"
    
    params = {
        "feature": "variation",
        "content-type": "application/json"
    }
    
    # Retry mechanism: attempt the request up to 'retries' times
    for attempt in range(retries):
        try:
            # Send the request to the Ensembl API
            response = requests.get(base_url, params=params, verify=False)
            
            # If the request is successful (status code 200), process the response
            if response.status_code == 200:
                data = response.json()
                # Extract SNP ID where all alleles are single characters (indicating SNV)
                snp_ids = [entry['id'] for entry in data if all(len(allele) == 1 for allele in entry.get('alleles', []))]
                # Return the first valid SNP ID if found; otherwise, return None
                return snp_ids[0] if snp_ids else None
            else:
                time.sleep(delay)  # Retry delay
        except requests.exceptions.RequestException as e:
            print(f"Network error on {chromosome}:{position} - {e}")
            time.sleep(delay)  # Retry delay if network error occurred

    # Return None if all retries failed
    return None

# Function to process each row in parallel
def process_row(row):
    chromosome = row['CHR'].replace('chr', '')  # Removing 'chr' prefix for API compatibility
    position = row['POS']
    gene_name = row['GNAME']
    sig_value = row['FEPVALUE']
    context_value = 'CD4_Mvalue09'  # As specified
    
    # Fetch SNP rs number
    snp_id = get_snp_rs_grch37(chromosome, position)
    
    # Return the processed data as a dictionary if SNP ID was retrieved
    if snp_id:
        return {
            'SNP': snp_id,
            'Chr': chromosome,
            'Pos': position,
            'Gene': gene_name,
            'Sig': sig_value,
            'Context': context_value
        }
    else:
        # Add to unrecognized SNPs if no ID was retrieved
        return {
            'SNP': None,
            'Chr': chromosome,
            'Pos': position,
            'Gene': gene_name,
            'Sig': sig_value,
            'Context': context_value
        }

# Read the BED file
file_path = '/Users/fvashegh/Library/CloudStorage/OneDrive-JNJ/Pi/eQTL/SynovialFibroblasts/eQTL/CD4_Mvalue09.bed'
cd4_eqtl_data = pd.read_csv(file_path, sep="\t") # nrows=3000
# Sort the DataFrame by 'FEPVALUE' in ascending order
cd4_eqtl_data = cd4_eqtl_data.sort_values(by='FEPVALUE', ascending=True)

# Initialize output lists
output_data = []
unrecognized_snps = []

# Start measuring time
start_time = time.time()

# Use ThreadPoolExecutor to process rows in parallel
with ThreadPoolExecutor(max_workers=100) as executor:
    # Create a future for each row
    futures = [executor.submit(process_row, row) for index, row in cd4_eqtl_data.iterrows()]
    
    # Collect results as they are completed
    for i, future in enumerate(as_completed(futures)):
        result = future.result()
        if result and result['SNP']:  # Only add recognized SNPs to output_data
            output_data.append(result)
        else:  # Add unrecognized SNPs to the unrecognized_snps list
            unrecognized_snps.append(result)
        
        # Print progress every 100 entries processed
        if (i + 1) % 1000 == 0:
            print(f"Processed {i + 1} entries...")

# Convert the recognized SNPs to a DataFrame
output_df = pd.DataFrame(output_data)

# Save the recognized SNPs as an RDS file for easy loading into R
output_file_path_rds = '/Users/fvashegh/Library/CloudStorage/OneDrive-JNJ/Pi/eQTL/SynovialFibroblasts/eQTL/processed_cd4_eqtl_data.RDS'
pyreadr.write_rds(output_file_path_rds, output_df)  # Save the DataFrame as RDS
print("Final SNP data saved to:", output_file_path_rds)

# End measuring time and calculate duration
end_time = time.time()
duration = end_time - start_time
print(f"Total runtime: {duration:.2f} seconds")










