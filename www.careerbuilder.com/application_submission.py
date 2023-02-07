import time
from selenium import webdriver

# input variables
job_text_field = "Software Engineer"
name = "John Doe"
phone = "555-555-1234"
email = "john.doe@email.com"
address = "1234 Main St"
city = "Anytown"
state = "CA"
zip_code = "12345"
resume_file = "/path/to/resume.pdf"

# list to store applied jobs
applied_jobs = []

# start a new instance of google chrome
driver = webdriver.Chrome()

# navigate to the careerbuilder.com homepage
driver.get("https://www.careerbuilder.com/")

# search for jobs with the keyword "Software Engineer"
search_field = driver.find_element_by_id("searchTerm")
search_field.clear()
search_field.send_keys(job_text_field)
search_button = driver.find_element_by_id("search-submit")
search_button.click()

# loop through each job listing
while True:
    try:
        # get the job title
        title = driver.find_element_by_xpath("//h2[@class='job-title']/a").text
        if title in applied_jobs:
            # move on to the next job if already applied
            next_button = driver.find_element_by_xpath("//a[@rel='next']")
            next_button.click()
            continue

        # click on the job title to view the job details
        driver.find_element_by_xpath("//h2[@class='job-title']/a").click()
        time.sleep(5)

        # click on the "Apply Now" button
        apply_button = driver.find_element_by_xpath("//a[@class='btn-apply']")
        apply_button.click()
        time.sleep(5)

        # fill in the form with personal information
        name_field = driver.find_element_by_xpath("//input[@name='fullName']")
        name_field.clear()
        name_field.send_keys(name)
        
        phone_field = driver.find_element_by_xpath("//input[@name='phone']")
        phone_field.clear()
        phone_field.send_keys(phone)
        
        email_field = driver.find_element_by_xpath("//input[@name='email']")
        email_field.clear()
        email_field.send_keys(email)

        # upload the resume
        resume_field = driver.find_element_by_xpath("//input[@type='file']")
        resume_field.send_keys(resume_file)

        # submit the form
        submit_button = driver.find_element_by_xpath("//button[@type='submit']")
        submit_button.click()
        time.sleep(5)
        
        # add this job to the list of applied jobs
        applied_jobs.append(title)

        # return to the
        driver.back()
        time.sleep(5)
    except Exception as e:
        # break the loop if an error occurs
        print(f"An error occurred: {e}")
        break
driver.quit()
