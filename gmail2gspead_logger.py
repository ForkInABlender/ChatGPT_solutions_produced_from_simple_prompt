import gspread
from oauth2client.service_account import ServiceAccountCredentials

# Authenticate to Google Drive
scope = ['https://spreadsheets.google.com/feeds',
         'https://www.googleapis.com/auth/drive']
credentials = ServiceAccountCredentials.from_json_keyfile_name('client_secret.json', scope)
gc = gspread.authorize(credentials)

# Open the Google Spreadsheet
spreadsheet = gc.open('Automated Email Log')
worksheet = spreadsheet.get_worksheet(0)

# Define the function to log new emails
def log_new_emails():
    result, data = mail.search(None, "ALL")
    spam_emails = data[0].split()
    for spam_email in spam_emails:
        result, data = mail.fetch(spam_email, "(RFC822)")
        raw_email = data[0][1].decode("utf-8")
        email_message = email.message_from_string(raw_email)
        sender = email_message["From"]
        subject = email_message["Subject"]
        message = email_message.get_payload()
        worksheet.append_row([sender, subject, message])

# Continuously check for new emails
while True:
    log_new_emails()
    time.sleep(30) # wait for 30 seconds before checking for new emails again


"""
Just logs content. Nothing special. Just adds a new colomn or row and adds  the data into it from your email.

Enjoy! 

"""
