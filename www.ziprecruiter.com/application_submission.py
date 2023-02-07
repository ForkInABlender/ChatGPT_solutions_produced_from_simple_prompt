from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import time

# initialize Chrome driver
driver = webdriver.Chrome()

# list of job titles to search for
job_list = ['Software Engineer', 'Data Scientist', 'Product Manager']

# loop through the job list
for job_title in job_list:
    # navigate to ZipRecruiter website
    driver.get("https://www.ziprecruiter.com/")

    # find the search bar and enter the job title
    search_bar = driver.find_element_by_id("home-search-field")
    search_bar.send_keys(job_title)

    # click the search button
    search_button = driver.find_element_by_id("search_submit_home")
    search_button.click()

    # wait for the results to load
    wait = WebDriverWait(driver, 10)
    jobs = wait.until(EC.presence_of_all_elements_located((By.CSS_SELECTOR, "a.job_link")))

    # loop through the list of jobs and apply to each one
    for job in jobs:
        job_link = job.get_attribute("href")
        driver.get(job_link)

        # check if application form exists
        try:
            apply_form = driver.find_element_by_css_selector("form#apply-form")
        except NoSuchElementException:
            print(f"No application form found for job: {job_title}")
            continue

        # fill in your personal information
        name_field = driver.find_element_by_css_selector("input#applicant_name")
        name_field.send_keys("Your Name")

        email_field = driver.find_element_by_css_selector("input#applicant_email")
        email_field.send_keys("your_email@example.com")

        phone_field = driver.find_element_by_css_selector("input#applicant_phone")
        phone_field.send_keys("1234567890")

        # attach your resume
        resume_field = driver.find_element_by_css_selector("input#resume")
        resume_field.send_keys("/path/to/your/resume.pdf")

        # submit the application form
        submit_button = driver.find_element_by_css_selector("input#submit_application")
        submit_button.click()

        # wait for the form to submit
        time.sleep(5)

# close the Chrome driver
driver.quit()
