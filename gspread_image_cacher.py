import gspread
from oauth2client.service_account import ServiceAccountCredentials
from PIL import Image

def image_to_cache(image_path, spreadsheet_title):
    scope = ['https://spreadsheets.google.com/feeds',
             'https://www.googleapis.com/auth/drive']
    credentials = ServiceAccountCredentials.from_json_keyfile_name('path/to/credentials.json', scope)
    client = gspread.authorize(credentials)

    image = Image.open(image_path)
    spreadsheet = client.create(spreadsheet_title)

    worksheet = spreadsheet.add_worksheet(title=image_path, rows=image.height, cols=image.width)

    for x in range(image.width):
        for y in range(image.height):
            pixel = image.getpixel((x, y))
            worksheet.update_cell(y+1, x+1, str(pixel))

def cache_to_image(spreadsheet_title, worksheet_title):
    scope = ['https://spreadsheets.google.com/feeds',
             'https://www.googleapis.com/auth/drive']
    credentials = ServiceAccountCredentials.from_json_keyfile_name('path/to/credentials.json', scope)
    client = gspread.authorize(credentials)

    spreadsheet = client.open(spreadsheet_title)
    worksheet = spreadsheet.worksheet(worksheet_title)

    image_height = worksheet.row_count
    image_width = worksheet.col_count

    image = Image.new('RGB', (image_width, image_height))

    for x in range(image_width):
        for y in range(image_height):
            cell_value = worksheet.cell(y+1, x+1).value
            image.putpixel((x, y), eval(cell_value))

    return image
