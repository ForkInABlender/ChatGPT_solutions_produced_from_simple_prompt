# import required modules
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import time

# setup input variables
job_text_field = "software engineer"
your_name = "John Doe"
phone_number = "1234567890"
email = "johndoe@email.com"
address = "123 Main St"
city = "Anytown"
state = "CA"
zip_code = "12345"

# setup the webdriver
options = webdriver.ChromeOptions()
options.add_argument('--disable-extensions')
options.add_argument('--profile-directory=Default')
options.add_argument("--incognito")
options.add_argument("--disable-plugins-discovery");
options.add_argument("--start-maximized")
driver = webdriver.Chrome(chrome_options=options)

# navigate to www.snagajob.com
driver.get("https://www.snagajob.com/")
time.sleep(5)

# search for jobs matching the `job_text_field`
search_box = driver.find_element_by_xpath("//input[@id='search']")
search_box.clear()
search_box.send_keys(job_text_field + Keys.RETURN)
time.sleep(5)

# loop through the job listings
while True:
    # find the job listings
    job_listings = driver.find_elements_by_xpath("//a[@class='job-title-link']")
    
    # break the loop if there are no more job listings
    if not job_listings:
        break
    
    # loop through the job listings
    for job in job_listings:
        # get the job title and link
        title = job.text
        link = job.get_attribute("href")
        
        # check if this job has already been applied for
        if title in applied_jobs:
            continue
        
        # navigate to the job listing page
        driver.get(link)
        time.sleep(5)
        
        # click the "Apply Now" button
        apply_now_button = driver.find_elements_by_xpath("//a[@class='apply-now-btn']")
        if apply_now_button:
            apply_now_button[0].click()
        else:
            continue
        
        # fill in the personal information form
        name_field = driver.find_element_by_xpath("//input[@name='full_name']")
        name_field.clear()
        name_field.send_keys(your_name)
        
        phone_field = driver.find_element_by_xpath("//input[@name='phone']")
        phone_field.clear()
        phone_field.send_keys(phone_number)
        
        email_field = driver.find_element_by_xpath("//input[@name='email']")
        email_field.clear()
        email_field.send_keys(email)
        
         # fill in the address information
        address_field = driver.find_element_by_xpath("//input[@name='address']")
        address_field.clear()
        address_field.send_keys(address)
        
        city_field = driver.find_element_by_xpath("//input[@name='city']")
        city_field.clear()
        city_field.send_keys(city)
        
        state_field = driver.find_element_by_xpath("//select[@name='state']")
        state_field.send_keys(state)
        
        zip_field = driver.find_element_by_xpath("//input[@name='zip']")
        zip_field.clear()
        zip_field.send_keys(zip_code)
        
        # submit the form
        submit_button = driver.find_element_by_xpath("//input[@value='Submit']")
        submit_button.click()
        time.sleep(5)
        
        # add this job to the list of applied jobs
        applied_jobs.append(title)

# close the webdriver
driver.quit()
