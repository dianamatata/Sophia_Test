
# libraries ------------
import time
import pandas as pd
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException
from selenium.webdriver.chrome.options import Options
from itertools import islice

# Set up Chrome options for headless mode
chrome_options = Options()
chrome_options.add_argument("--headless")  # Run in headless mode
chrome_options.add_argument("--no-sandbox")  # Optional: avoid some sandboxing issues in certain environments
chrome_options.add_argument("--disable-dev-shm-usage")  # Optional: solve issues in Docker environments



# functions ------------
def extract_data_from_website(driver):
    """
    Extracts the gene name and table data from the website and returns them.
    The data is extracted using Selenium and returned as a tuple (gene_name, dataframe).
    """
    wait = WebDriverWait(driver, 10)

    # Wait until the gene name element is visible
    gene_element = wait.until(EC.visibility_of_element_located((By.ID, "gene_symbol")))
    gene_name = gene_element.text
    # gene_element = driver.find_element(By.ID, "gene_symbol")
    # gene_name = gene_element.find_element(By.ID, "gene_symbol").text
    print(f"Gene Name: {gene_name}")

    # Extract all rows from the table and store them in a DataFrame
    table_rows = driver.find_elements(By.XPATH, "//tr")
    table_data = []

    # Collect table data
    for row in table_rows:
        cells = row.find_elements(By.TAG_NAME, "td")
        row_data = [cell.text for cell in cells]
        if row_data:  # Ensure the row has content
            table_data.append(row_data)

    # Convert the table data into a DataFrame
    df = pd.DataFrame(table_data)
    df = df[df[0].notna() & (df[0] != '')]
    df = df.drop(columns=[2, 3, 4])
    tdf = df.T
    tdf.columns = tdf.iloc[0]
    tdf = tdf.drop(tdf.index[0])

    return gene_name, tdf

# File paths ------------------
formatted_entries_file = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/formatted_vcf_entries.txt"
output_file = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/variant_NM_VCF.csv"
# formatted_entries_file = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/test_entries.txt"
failed_variants_log = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/failed_variants_log.txt"

# main ------------------

# Track the last processed line
last_processed_line_file = '/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/last_processed_line.txt'
output_mobidetails = '/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/output/output_mobidetails'

# Check if the file exists, and if so, read the last processed line
try:
    with open(last_processed_line_file, 'r') as f:
        last_processed_line = int(f.read().strip())
except FileNotFoundError:
    last_processed_line = 0  # If the file does not exist, start from the beginning


# Initialize the webdriver
driver = webdriver.Chrome(options=chrome_options)
driver.get("https://mobidetails.iurc.montp.inserm.fr/MD/api/variant/5/browser/")

# DEBUG -----
# Read the file and get the variant from line 27
with open(formatted_entries_file, 'r') as f:
    lines = f.readlines()  # Read all lines
    variant = lines[26].strip()  # Line 27 is index 26 (since indexing starts from 0)
    print(f"Variant from line 27: {variant}")
# DEBUG -----

# Initialize the iteration counter and tdf1
tdf1 = pd.DataFrame()

last_processed_line = line_num -1
iteration = 7

# Read the formatted VCF entries from the file
with open(formatted_entries_file, 'r') as f:
    # Skip the lines that have already been processed
    for _ in range(last_processed_line):
        next(f)

    # Loop through each line (variant) in the file
    for line_num, line in enumerate(f, start=last_processed_line + 1):
        variant = line.strip()  # Get the current variant
        print(f"Processing variant {line_num}: {variant}")

        # Locate the input field and fill it with the value
        input_element = WebDriverWait(driver, 30).until(
            EC.presence_of_element_located((By.ID, "engine")))
        input_element.clear()
        input_element.send_keys(variant)
        search_button = driver.find_element(By.ID, "submit_a")
        search_button.click()

        # maybe we have to pick variant in several transcripts
        case_1_xpath = "//td[@class='w3-left-align dtr-control' and text()='HGNC gene symbol (ID):']"
        case_2_xpath = "//li[.//em[starts-with(text(), 'NM_')]]"
        try:
            # Try to find the first case
            WebDriverWait(driver, 5).until(EC.presence_of_element_located((By.XPATH, case_1_xpath)))
        except TimeoutException:
            # If case 1 fails, try case 2
            try:
                WebDriverWait(driver, 5).until(EC.presence_of_element_located((By.XPATH, case_2_xpath)))
                # Extract the full variant name from the <em> tag inside the <li> element
                variant_element = driver.find_element(By.XPATH, case_2_xpath)
                # Extract the URL from the onclick attribute
                onclick_script = variant_element.get_attribute("onclick")
                # Combine base URL with the relative URL
                full_url = "https://mobidetails.iurc.montp.inserm.fr" + onclick_script.split("window.open('")[1].split("')")[0]
                # Navigate to the full URL directly
                driver.get(full_url)
            except TimeoutException:
                # Log the failed variant
                with open(failed_variants_log, 'a') as log_file:
                    log_file.write(f"Failed to process variant {variant} on line {line_num}\n")
                print(f"Variant {variant} failed on line {line_num}. Skipping to next.")
                continue  # Skip to the next variant

        time.sleep(3)     # Wait for a moment before the next search to avoid overloading the server

        # RUN FUNCTION TO GET DATA
        gene_name, df_var = extract_data_from_website(driver)
        tdf1 = pd.concat([tdf1, df_var], axis=0, ignore_index=True, sort=False)

        if line_num % 50 == 0:  # Check if the number of rows is a multiple of 50, to not write all the time
            output_filename = f"{output_mobidetails}_{iteration}.txt"
            tdf1.to_csv(output_filename, index=False)
            print(f"Saved to {output_filename}")
            tdf1 = pd.DataFrame() # reinitialize
            iteration += 1
            last_processed_line = line_num
            with open(last_processed_line_file, 'w') as f:
                f.write(str(last_processed_line))


# Close the driver
# driver.quit()

# check
output_filename = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/output_mobidetails_2.txt"
df = pd.read_csv(output_filename)
print(df.shape)
# Check the contents of the DataFrame
print(df.head())

# Close the browser
driver.quit()
