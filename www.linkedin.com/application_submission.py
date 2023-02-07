from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import time

# list of jobs to apply for
jobs = ["Software Engineer", "Data Scientist", "DevOps Engineer", "Product Manager"]

# personal information to fill in
name = "John Doe"
phone = "123-456-7890"
email = "johndoe@example.com"
address = "1234 Main St, Anytown USA 12345"

# initialize the webdriver
driver = webdriver.Chrome()

# loop through each job in the list
for job in jobs:
    # navigate to LinkedIn job search page
    driver.get("https://www.linkedin.com/jobs/")
    
    # wait for the page to load
    time.sleep(2)
    
    # search for the job
    search_field = driver.find_element_by_xpath("//input[@role='combobox']")
    search_field.send_keys(job + Keys.RETURN)
    
    # wait for the search results to load
    time.sleep(2)
    
    # click the first job result
    first_job = driver.find_element_by_xpath("//a[@data-link-to='job']")
    first_job.click()
    
    # wait for the job page to load
    time.sleep(2)
    
    # click the "Apply" button
    apply_button = driver.find_element_by_xpath("//button[text()='Apply']")
    apply_button.click()
    
    # wait for the application form to load
    time.sleep(2)
    
    # fill in the personal information
    name_field = driver.find_element_by_xpath("//input[@id='name']")
    name_field.send_keys(name)
    
    phone_field = driver.find_element_by_xpath("//input[@id='phone']")
    phone_field.send_keys(phone)
    
    email_field = driver.find_element_by_xpath("//input[@id='email']")
    email_field.send_keys(email)
    
    address_field = driver.find_element_by_xpath("//textarea[@id='address']")
    address_field.send_keys(address)
    
    # submit the application
    submit_button = driver.find_element_by_xpath("//button[text()='Submit']")
    submit_button.click()
    
    # wait for the confirmation page to load
    time.sleep(2)
    
    # navigate back to the job search page
    driver.get("https://www.linkedin.com/jobs/")

# close the browser
driver.quit()
