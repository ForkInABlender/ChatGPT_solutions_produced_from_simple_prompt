from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By

# setup chrome driver
driver = webdriver.Chrome()

# get the applicant's information
name = input("What is your name? ")
phone_number = input("What is your phone number? ")
email = input("What is your email address? ")
address = input("What is your address? ")

# navigate to Indeed's website
driver.get("https://www.indeed.com/")

# find the "Job Title" field and fill it in
job_title_field = driver.find_element_by_id("text-input-what")
job_title_field.send_keys("Software Engineer")

# find the "Location" field and fill it in
location_field = driver.find_element_by_id("text-input-where")
location_field.send_keys("San Francisco, CA")

# find the "Find Jobs" button and click it
find_jobs_button = driver.find_element_by_xpath("//button[text()='Find Jobs']")
find_jobs_button.click()

while True:
    # wait for the job listings to load
    job_listings = WebDriverWait(driver, 10).until(EC.presence_of_all_elements_located((By.XPATH, "//div[@data-tn-component='organicJob']")))

    # loop through each job listing
    for job_listing in job_listings:
        # extract the job title
        job_title = job_listing.find_element_by_xpath(".//h2").text

        # only apply for jobs that match the job_title_field string
        if job_title != job_title_field:
            continue

        # check if you've already applied for this job
        applied_label = job_listing.find_elements_by_xpath(".//span[text()='Applied']")
        if applied_label:
            continue

        # click the job listing
        job_listing.click()

        # find the "Apply" button and click it
        apply_button = driver.find_element_by_xpath("//button[text()='Apply']")
        apply_button.click()

        # fill in the applicant's information
        name_field = driver.find_element_by_name("fullName")
        name_field.send_keys(name)

        phone_number_field = driver.find_element_by_name("phone")
        phone_number_field.send_keys(phone_number)

        email_field = driver.find_element_by_name("email")
        email_field.send_keys(email)

        address_field = driver.find_element_by_name("address")
        address_field.send_keys(address)

        # submit the application
        submit_button = driver.find_element_by_xpath("//button[text()='Submit']")
                submit_button.click()
        
        # mark the job as applied
        job_listing.find_element_by_xpath(".//span[text()='Applied']").click()
        
        # go back to the job search results page
        driver.back()

# close the browser window
driver.quit()
