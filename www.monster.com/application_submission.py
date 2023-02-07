from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import time

# Create a list of job openings you want to apply to
job_list = ["Software Engineer", "Data Scientist", "Product Manager"]

# Information for the job application
name = "John Doe"
email = "johndoe@email.com"
phone = "555-555-5555"
resume = "/path/to/resume.pdf"

# Start Chrome webdriver
driver = webdriver.Chrome()

# Loop through the list of jobs
for job in job_list:
    # Navigate to Monster.com
    driver.get("https://www.monster.com/")

    # Find the search bar and input the desired job title
    search_bar = driver.find_element_by_name("q")
    search_bar.send_keys(job)
    search_bar.send_keys(Keys.RETURN)

    # Wait for search results to load
    time.sleep(5)

    # Find the first job result and click on it
    first_result = driver.find_element_by_css_selector("#ResultsContainer .js_result_container:first-of-type .jobTitle")
    first_result.click()

    # Wait for job description to load
    time.sleep(5)

    # Check if you have already applied for this job
    applied = False
    try:
        driver.find_element_by_css_selector(".btn.btn-secondary.already-applied")
        applied = True
    except:
        pass

    # If you haven't applied, continue with the application process
    if not applied:
        # Find the "Apply" button and click on it
        apply_button = driver.find_element_by_css_selector(".btn.btn-primary.js-btn-apply")
        apply_button.click()

        # Wait for the application form to load
        time.sleep(5)

        # Fill out the form with your information
        name_field = driver.find_element_by_name("FullName")
        name_field.send_keys(name)

        email_field = driver.find_element_by_name("Email")
        email_field.send_keys(email)

        phone_field = driver.find_element_by_name("Phone")
        phone_field.send_keys(phone)

        resume_field = driver.find_element_by_name("Resume")
        resume_field.send_keys(resume)

        # Submit the application
        submit_button = driver.find_element_by_css_selector(".btn.btn-primary.js-btn-submit")
        submit_button.click()

        # Wait for confirmation page to load
        time.sleep(5)

        # Close the application confirmation page
        driver.back()

# Quit the webdriver
driver.quit()
