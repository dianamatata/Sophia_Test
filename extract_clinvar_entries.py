from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import pandas as pd

# Initialize WebDriver
driver = webdriver.Chrome()  # Update with the path to your ChromeDriver if not in PATH
url = "https://www.ncbi.nlm.nih.gov/clinvar/variation/100094/?oq=%22NM_000110.4(DPYD):c.1601G%3EA%20(p.Ser534Asn)%22%5BVARNAME%5D&m=NM_000110.4(DPYD):c.1601G%3EA%20(p.Ser534Asn)%3Fterm=NM_000110.4(DPYD):c.1601G%3EA%20(p.Ser534Asn)"
driver.get(url)

try:
    # Wait for the table to load
    table = WebDriverWait(driver, 10).until(
        EC.presence_of_element_located((By.CSS_SELECTOR,
                                        "table.stickyheaders.usa-table-borderless.table-top-border.align-top.conditions-germline-list.expandable"))
    )

    # Extract rows of the table
    rows = table.find_elements(By.TAG_NAME, "tr")

    # Initialize a list to store extracted data
    data = []
    for row in rows:
        cells = row.find_elements(By.TAG_NAME, "td")
        row_data = [cell.text for cell in cells]
        if row_data:  # Skip rows without data
            data.append(row_data)

    # Create DataFrame
    headers = [th.text.strip() for th in table.find_elements(By.TAG_NAME, "th")]
    df = pd.DataFrame(data, columns=headers)

    # Print DataFrame
    print(df)

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

    # Output the dictionary
    print(classification_dict)


# TODO we want a function that iterates over all the variants that have clinvar LP or P , and add it to the mobidetails df column



# search over var

    search_field = driver.find_element(By.ID, "search-field")
    search_field.send_keys("eee")
    submit_button = driver.find_element(By.XPATH, "//button[@type='submit']")
    submit_button.click()



    driver.quit()
