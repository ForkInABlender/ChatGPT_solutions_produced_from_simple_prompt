from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
import time

# create instance of chrome webdriver
driver = webdriver.Chrome()

# list of jobs to apply for
jobs = ["Software Engineer", "Data Scientist", "Product Manager", "DevOps Engineer"]

for job in jobs:
    # navigate to glassdoor website
    driver.get("https://www.glassdoor.com/")

    # wait for the page to load
    time.sleep(3)

    # find the search bar element and enter job_text_field string as the search keyword
    search_bar = driver.find_element(By.ID, "sc.keyword")
    search_bar.send_keys(job)
    search_bar.send_keys(Keys.RETURN)

    # wait for the search results to load
    time.sleep(3)

    # find the first job in the search results and click on it
    first_job = driver.find_element(By.XPATH, "//li[1]//div[1]//div[1]//a[1]")
    first_job.click()

    # wait for the job page to load
    time.sleep(3)

    # find the apply button and click on it
    apply_button = driver.find_element(By.XPATH, "//button[@data-test='job-apply']")
    apply_button.click()

    # wait for the application form to load
    time.sleep(3)

    # fill in personal information
    name_field = driver.find_element(By.ID, "firstName")
    name_field.send_keys(full_name)

    phone_field = driver.find_element(By.ID, "phone")
    phone_field.send_keys(phone_number)

    email_field = driver.find_element(By.ID, "email")
    email_field.send_keys(email)

    address_field = driver.find_element(By.ID, "address")
    address_field.send_keys(address)

    # upload resume
    upload_resume_button = driver.find_element(By.ID, "file-upload")
    upload_resume_button.send_keys(resume_file_path)

    # submit the application
    submit_button = driver.find_element(By.XPATH, "//button[@type='submit']")
    submit_button.click()

# close the browser
driver.quit()
