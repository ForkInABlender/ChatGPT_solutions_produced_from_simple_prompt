from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import time

# specify the path to the Chrome driver
driver = webdriver.Chrome(executable_path=r"/path/to/chromedriver")

# open USAJobs.gov in the browser
driver.get("https://www.usajobs.gov/")

# wait for the page to load
time.sleep(5)

# find the search input field and enter the job text field
search_input = driver.find_element_by_id("search-input-field")
job_text_field = "software developer"
search_input.send_keys(job_text_field)

# find the search button and click it
search_button = driver.find_element_by_id("search-submit-button")
search_button.click()

# wait for the search results page to load
time.sleep(5)

# get a list of all job listings on the page
job_listings = driver.find_elements_by_xpath("//ul[@id='searchResultsList']/li")

# loop through each job listing
for job_listing in job_listings:
    # get the title and location of the job
    title = job_listing.find_element_by_xpath(".//h2").text
    location = job_listing.find_element_by_xpath(".//div[@class='job-location']").text
    
    # check if the job title and location match the job text field
    if job_text_field in title and job_text_field in location:
        # open the job listing
        job_listing.find_element_by_xpath(".//a[@class='job-title-link']").click()
        
        # wait for the job listing page to load
        time.sleep(5)
        
        # check if the Apply button is present
        apply_button = driver.find_elements_by_xpath("//button[text()='Apply']")
        if apply_button:
            # fill in the personal information form
            driver.find_element_by_id("Resume").send_keys("/path/to/resume.pdf")
            driver.find_element_by_id("firstName").send_keys("John")
            driver.find_element_by_id("lastName").send_keys("Doe")
            driver.find_element_by_id("contactInfo.email").send_keys("john.doe@example.com")
            driver.find_element_by_id("contactInfo.phone").send_keys("1234567890")
            driver.find_element_by_id("contactInfo.address.addressLine1").send_keys("123 Main St")
            driver.find_element_by_id("contactInfo.address.city").send_keys("Washington")
            driver.find_element_by_id("contactInfo.address.state").send_keys("DC")
            driver.find_element_by_id("contactInfo.address.zipCode").send_keys("20000")
            # submit the form
            apply_button[0].click()
            
            # wait for the application confirmation page to load
            time.sleep(5)
            
            # check if the application was successfully submitted
            confirmation_text = driver.find_element_by_xpath("//div[@class='usa-alert-body']").text
            if "Your application has been submitted" in confirmation_text:
                print("Successfully applied for job:", title)
            else:
                print("Failed to apply for job:", title)
                
            # go back to the job search results page
            driver.back()
            time.sleep(5)

# close the browser
driver.quit()
