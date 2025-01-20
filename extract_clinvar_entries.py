from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import pandas as pd
pd.set_option('display.max_columns', None)

# apply code to subset_clinvar_pathogenic
subset_clinvar_pathogenic = mobi_data[
    mobi_data['simplified_clinvar'].isin([
        'drug_response',
        'Likely_pathogenic',
        'Pathogenic/Likely_pathogenic',
        'Pathogenic'
    ])
]

print(subset_clinvar_pathogenic.shape)
print(subset_clinvar_pathogenic['simplified_clinvar'])



# Initialize WebDriver
driver = webdriver.Chrome()  # Update with the path to your ChromeDriver if not in PATH
url = "https://www.ncbi.nlm.nih.gov/clinvar/variation/100094/?oq=%22NM_000110.4(DPYD):c.1601G%3EA%20(p.Ser534Asn)%22%5BVARNAME%5D&m=NM_000110.4(DPYD):c.1601G%3EA%20(p.Ser534Asn)%3Fterm=NM_000110.4(DPYD):c.1601G%3EA%20(p.Ser534Asn)"
driver.get(url)

# Loop through subset_clinvar_pathogenic
# row_df = subset_clinvar_pathogenic.iloc[2]
for index, row_df in subset_clinvar_pathogenic.iterrows():
    print(row_df)
    value = row_df.get('dbSNP rsid:', None)
    # Locate the search field and input the search term
    search_field = driver.find_element(By.ID, "search-field")
    search_field.clear()  # Clear any previous text
    search_field.send_keys(value)

    # Locate and click the submit button
    submit_button = driver.find_element(By.XPATH, "//button[@type='submit']")
    submit_button.click()

    # get_classification_dict(driver, subset_clinvar_pathogenic)
    # Wait for the table to load
    table = WebDriverWait(driver, 10).until(
        EC.presence_of_element_located((By.CSS_SELECTOR,
                                        "table.stickyheaders.usa-table-borderless.table-top-border.align-top.conditions-germline-list.expandable"))
    )
    rows = table.find_elements(By.TAG_NAME, "tr")     # Extract rows of the table

    data = []
    for row in rows:
        cells = row.find_elements(By.TAG_NAME, "td")
        row_data = [cell.text for cell in cells]
        if row_data:  # Skip rows without data
            data.append(row_data)

    headers = [th.text.strip() for th in table.find_elements(By.TAG_NAME, "th")]
    df = pd.DataFrame(data, columns=headers)
    # Keep only the required columns and rename them
    df = df.iloc[:, [0, 1]].rename(columns={df.columns[0]: "Condition", df.columns[1]: "Classification"})

    classification_dict = {}
    # Loop through the 'Classification' column to process the values
    # Iterate over each entry in the 'Classification' column
    for entry in df['Classification']:
        # Split the entry into name and value based on '('
        parts = entry.split('(')
        name = parts[0].strip()  # Extract the name and remove extra whitespace
        count = int(parts[1].strip(')')) if len(parts) > 1 else 0  # Extract the value or set to 0 if not present
        # Increment the dictionary value or add the key if it doesn't exist
        classification_dict[name] = classification_dict.get(name, 0) + count

    print(classification_dict)
    classification_dict_str = ', '.join([f"{key}: {value}" for key, value in classification_dict.items()])

    subset_clinvar_pathogenic.loc[index, 'Clinvar_entries'] = classification_dict_str

driver.quit()


# subset_clinvar_pathogenic.rename(columns={'CliVar entries': 'Clinvar_entries'}, inplace=True)

subset_clinvar_pathogenic.to_csv('/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/subset_clinvar_pathogenic.txt', index=True, sep='\t')
subset_clinvar_pathogenic = pd.read_csv(
    '/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/subset_clinvar_pathogenic.txt',
    sep='\t',
    index_col=0  # This makes the first column the index
)

# incorporate it in mobi_data
subset_clinvar_pathogenic['Clinvar_entries']
mobi_data['Clinvar_entries'] = ""
# Step 2: Incorporate the processed results back into mobi_data
mobi_data.loc[subset_clinvar_pathogenic.index, 'Clinvar_entries'] = subset_clinvar_pathogenic['Clinvar_entries']

#  we want a function that iterates over all the variants that have clinvar LP or P , and add it to the mobidetails df column


print(mobi_data.shape) # (2406, 69)
print(mobi_data['Clinvar_entries'].value_counts())

# save data -----------
mobi_data.to_csv('/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/mobi_data_omim_splice_clinvarentries.txt', index=False, sep='\t')
mobi_data = pd.read_csv("/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/mobi_data_omim_splice_clinvarentries.txt", sep='\t')


