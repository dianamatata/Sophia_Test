import pandas as pd
from selenium import webdriver
from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import time

# Set up Chrome options
options = Options()
options.headless = False  # Set to True if you don't want the browser window to open

# Set up WebDriver
driver = webdriver.Chrome(service=Service(ChromeDriverManager().install()), options=options)

# Read HPO IDs from the TSV file
hpo_ids = []
with open('/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/patient_phenotype.tsv', 'r') as f:
    hpo_ids = [line.strip() for line in f.readlines()]

# Initialize a list to store the results
results = []
hpo_id = 'HP:0000365'
# Iterate over each HPO ID, open the URL and extract the term name and gene associations
for hpo_id in hpo_ids:
    url = f'https://hpo.jax.org/browse/term/{hpo_id}'
    driver.get(url)

    # Wait for the page to load (adjust if necessary)
    driver.implicitly_wait(5)  # Wait for 5 seconds
    try:
        # Extract the term name
        term_name = driver.find_element(By.CSS_SELECTOR, '.item-title').text
    except Exception as e:
        term_name = 'Not Found'  # If the term is not found, mark as 'Not Found'

    # Click the "Gene Associations [Inferred]" tab
    try:
        gene_associations_button = WebDriverWait(driver, 10).until(
            EC.element_to_be_clickable(
                (By.XPATH, '//div[@class="mat-tab-label-content" and text()="Gene Associations [Inferred]"]'))
        )
        gene_associations_button.click()
    except Exception as e:
        print(f"Gene Associations button not found for HPO ID {hpo_id}")

    # Wait for the table to load
    WebDriverWait(driver, 10).until(
        EC.presence_of_element_located(
            (By.CSS_SELECTOR, '#mat-tab-content-0-1 .associations.gene-association mat-table'))
    )

    # Extract the gene association table
    try:
        table = driver.find_element(By.CSS_SELECTOR, '#mat-tab-content-0-1 .associations.gene-association mat-table')
        rows = table.find_elements(By.TAG_NAME, 'mat-row')

        # Extract data from the rows
        table_data = []
        for row in rows:
            cells = row.find_elements(By.TAG_NAME, 'mat-cell')
            table_data.append(cells[1].text)
    except Exception as e:
        table_data = []

    # Store the HPO ID, term name, and associated gene data
    results.append({'HPO_ID': hpo_id, 'Term_Name': term_name, 'Gene_Associations': table_data})

    # Optional: Sleep to avoid overloading the server
    time.sleep(1)

df = pd.DataFrame(results)
df.to_csv('/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/hpo_terms_and_gene_associations.csv', index=False)


# or save in json file
import json

# Save to a JSON file
with open('/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/hpo_terms_and_gene_associations.json', 'w') as json_file:
    json.dump(results, json_file, indent=4)

# Load the data later
with open('/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/hpo_terms_and_gene_associations.json', 'r') as json_file:
    loaded_data = json.load(json_file)


# Print the results
print(df)

# Close the driver
driver.quit()

df[['HPO_ID', 'Term_Name']].to_csv('/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/hpo_summary.txt', index=False, sep='\t')
