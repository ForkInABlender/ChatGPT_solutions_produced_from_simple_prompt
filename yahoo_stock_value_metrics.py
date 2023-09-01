#Dylan Kenneth Eliot & GPT-4-plugins

"""
This is a class that goes through the values on yahoo finances. Currently, the example only monitors 1 value on yahoo finances and compares it to the last value it had previously stored in memory. This also can be setup to watch multiple
 stock value metrics. Now, this is a template of what you can do. Most of it was simple to get working. If you need those metrics, uncomment or add in by html tag name and the thing you need it to identify as seen in similar examples within this
  class. This can also be adapted to other websites also listing stock price values. 









"""

import requests
from bs4 import BeautifulSoup

class yahoo_Stock_metrics_class:
	def __init__(self, url):
		self.headers= {
			"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
		}
		self.url=url
		#
	#
	def safe_extract(self, soup_obj, tag, class_name, alternative_class_name=None):
		result = soup_obj.find_all(tag, class_=class_name)
		if not result and alternative_class_name:
			result = soup_obj.find_all(tag, class_=alternative_class_name)
		return result[0].text if result else None

	def get_metrics(self):
		response = requests.get(self.url, headers=self.headers)
		soup = BeautifulSoup(response.content, 'html.parser')
		current_price = self.safe_extract(soup, "fin-streamer", "Fw(b) Fz(36px) Mb(-4px) D(ib)")
		# Extracting the data using the safe_extract function
		#current_price = safe_extract(soup, "fin-streamer", "Fw(b) Fz(36px) Mb(-4px) D(ib)")
		#change = safe_extract(soup, "fin-streamer", "Fw(500) Pstart(8px) Fz(24px)")
		#previous_close = soup.find(string="Previous Close").find_next("td").text if soup.find(string="Previous Close") else None
		#_open = soup.find(string="Open").find_next("td").text if soup.find(string="Open") else None
		#day_range = soup.find(string="Day's Range").find_next("td").text if soup.find(string="Day's Range") else None
		#week_range = soup.find(string="52 Week Range").find_next("td").text if soup.find(string="52 Week Range") else None
		#start_date = soup.find(string="Start Date").find_next("td").text if soup.find(string="Start Date") else None
		#market_cap = soup.find(string="Market Cap").find_next("td").text if soup.find(string="Market Cap") else None
		#circulating_supply = soup.find(string="Circulating Supply").find_next("td").text if soup.find(string="Circulating Supply") else None
		#volume = soup.find(string="Volume").find_next("td").text if soup.find(string="Volume") else None
		#volume_24hr = soup.find(string="Volume (24hr)").find_next("td").text if soup.find(string="Volume (24hr)") else None
		metrics = {
  	'Current Price': current_price#,
    #'Change': change,
    #'Previous Close': previous_close,
    #'Open': _open,
    #'Day\'s Range': day_range,
    #'52 Week Range': week_range,
    #'Start Date': start_date,
    #'Market Cap': market_cap,
    #'Circulating Supply': circulating_supply,
    #'Volume': volume,
    #'Volume (24hr)': volume_24hr
		}
		return metrics



urls = ["https://finance.yahoo.com/quote/BTC-USD?p=BTC-USD"]*100000

previous_metrics = None
for url in urls:
    current_metrics = yahoo_Stock_metrics_class(url).get_metrics()
    if previous_metrics and current_metrics['Current Price'] != previous_metrics['Current Price']:
        print(f"Change detected in stock value from {previous_metrics['Current Price']} to {current_metrics['Current Price']}")
    previous_metrics = current_metrics
