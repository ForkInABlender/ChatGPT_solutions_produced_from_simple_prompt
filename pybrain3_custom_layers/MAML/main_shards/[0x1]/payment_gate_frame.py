# Dylan Kenneth Eliot

"""
Basic payment gateway via flask & stripe


"""


from flask import Flask, request, render_template_string, redirect, url_for
import stripe

app = Flask(__name__)

# Stripe API configuration
stripe.api_key = ''

html_form = '''
<!DOCTYPE html>
<html lang="en">
<head>
		<meta charset="UTF-8">
		<meta name="viewport" content="width=device-width, initial-scale=1.0">
		<title>Option Selection Form</title>
		<script type="text/javascript" src="https://cdn.jsdelivr.net/npm/brython@3.9.5/brython.min.js"></script>
</head>
<body onload="brython()">
		<form action="/create-checkout-session" method="post" onsubmit="return validate_selection()">
				<label for="subscription_pack">Choose Subscription Pack (1-6):</label><br>
				<select id="subscription_pack" name="subscription_pack">
						<option value="1">Subscription Pack 1</option>
						<option value="2">Subscription Pack 2</option>
						<option value="3">Subscription Pack 3</option>
						<option value="4">Subscription Pack 4</option>
						<option value="5">Subscription Pack 5</option>
						<option value="6">Subscription Pack 6</option>
				</select><br><br>

				<label for="options">Choose one or more options:</label><br>
				<input type="checkbox" id="option1" name="option" value="option1" onchange="calculate_total()">
				<label for="option1"> Combination 1: 6 × C1 - $3.42</label><br>
				<input type="checkbox" id="option2" name="option" value="option2" onchange="calculate_total()">
				<label for="option2"> Combination 2: 1 × C1 + 1 × C2 - $3.83</label><br>
				<input type="checkbox" id="option3" name="option" value="option3" onchange="calculate_total()">
				<label for="option3"> Combination 3: 2 × C1 + 1 × C2 - $4.40</label><br>
				<input type="checkbox" id="option4" name="option" value="option4" onchange="calculate_total()">
				<label for="option4"> Combination 4: 3 × C1 + 1 × C2 - $4.97</label><br>
				<input type="checkbox" id="option5" name="option" value="option5" onchange="calculate_total()">
				<label for="option5"> Combination 5: 4 × C1 + 1 × C2 - $5.54</label><br>
				<input type="checkbox" id="option6" name="option" value="option6" onchange="calculate_total()">
				<label for="option6"> Combination 6: 5 × C1 + 1 × C2 - $6.11</label><br>
				<input type="checkbox" id="option7" name="option" value="option7" onchange="calculate_total()">
				<label for="option7"> Combination 7: 1 × C1 + 1 × C3 - $6.52</label><br>
				<input type="checkbox" id="option8" name="option" value="option8" onchange="calculate_total()">
				<label for="option8"> Combination 8: 2 × C1 + 1 × C3 - $7.09</label><br>
				<input type="checkbox" id="option9" name="option" value="option9" onchange="calculate_total()">
				<label for="option9"> Combination 9: 2 × C1 + 2 × C2 - $7.66</label><br>
				<input type="checkbox" id="option10" name="option" value="option10" onchange="calculate_total()">
				<label for="option10"> Combination 10: 3 × C1 + 1 × C3 - $7.66</label><br>
				<input type="checkbox" id="option11" name="option" value="option11" onchange="calculate_total()">
				<label for="option11"> Combination 11: 4 × C1 + 1 × C3 - $8.23</label><br>
				<input type="checkbox" id="option12" name="option" value="option12" onchange="calculate_total()">
				<label for="option12"> Combination 12: 5 × C1 + 1 × C3 - $8.80</label><br>
				<input type="checkbox" id="option13" name="option" value="option13" onchange="calculate_total()">
				<label for="option13"> Combination 13: 3 × C1 + 1 × C4 - $10.35</label><br>
				<input type="checkbox" id="option14" name="option" value="option14" onchange="calculate_total()">
				<label for="option14"> Combination 14: 1 × C1 + 2 × C3 - $12.47</label><br>
				<input type="checkbox" id="option15" name="option" value="option15" onchange="calculate_total()">
				<label for="option15"> Combination 15: 4 × C2 - $13.04</label><br>
				<input type="checkbox" id="option16" name="option" value="option16" onchange="calculate_total()">
				<label for="option16"> Combination 16: 3 × C1 + 1 × C5 - $13.04</label><br>
				<input type="checkbox" id="option17" name="option" value="option17" onchange="calculate_total()">
				<label for="option17"> Combination 17: 4 × C1 + 1 × C4 - $10.92</label><br>
				<input type="checkbox" id="option18" name="option" value="option18" onchange="calculate_total()">
				<label for="option18"> Combination 18: 5 × C1 + 1 × C4 - $11.49</label><br>
				<input type="checkbox" id="option19" name="option" value="option19" onchange="calculate_total()">
				<label for="option19"> Combination 19: 5 × C1 + 1 × C5 - $14.18</label><br>
				<input type="checkbox" id="option20" name="option" value="option20" onchange="calculate_total()">
				<label for="option20"> Combination 20: 1 × C3 + 1 × C5 + 1 × C1 - $17.85</label><br>
				<input type="checkbox" id="option21" name="option" value="option21" onchange="calculate_total()">
				<label for="option21"> Combination 21: 2 × C1 + 1 × C7 - $17.85</label><br>
				<input type="checkbox" id="option22" name="option" value="option22" onchange="calculate_total()">
				<label for="option22"> Combination 22: 5 × C1 + 1 × C6 - $16.87</label><br>
				<input type="checkbox" id="option23" name="option" value="option23" onchange="calculate_total()">
				<label for="option23"> Combination 23: 5 × C1 + 1 × C7 - $19.56</label><br>
				<input type="checkbox" id="option24" name="option" value="option24" onchange="calculate_total()">
				<label for="option24"> Combination 24: 1 × C1 + 1 × C8 - $19.99</label><br><br>
				<p id="total_cost">Total Cost: $0.00</p>
				<input type="submit" value="Submit">
		</form>

		<script type="text/python">
		from browser import document, alert

		option_costs = {
				"option1": 3.42,
				"option2": 3.83,
				"option3": 4.40,
				"option4": 4.97,
				"option5": 5.54,
				"option6": 6.11,
				"option7": 6.52,
				"option8": 7.09,
				"option9": 7.66,
				"option10": 7.66,
				"option11": 8.23,
				"option12": 8.80,
				"option13": 10.35,
				"option14": 12.47,
				"option15": 13.04,
				"option16": 13.04,
				"option17": 10.92,
				"option18": 11.49,
				"option19": 14.18,
				"option20": 17.85,
				"option21": 17.85,
				"option22": 16.87,
				"option23": 19.56,
				"option24": 19.99
		}

		def calculate_total():
				total_cost = 0
				for option in document.select("input[type=checkbox]"):
						if option.checked:
								total_cost += option_costs[option.value]
				document["total_cost"].text = f"Total Cost: ${total_cost:.2f}"

				# Alert if total exceeds $20
				if total_cost > 20:
						alert(f"Warning: Selection exceeds $20. Total cost: ${total_cost:.2f}.")

		def validate_selection():
				total_cost = 0
				for option in document.select("input[type=checkbox]"):
						if option.checked:
								total_cost += option_costs[option.value]

				if total_cost > 20:
						alert(f"Selection exceeds $20. Total cost: ${total_cost:.2f}. Please choose fewer options.")
						return False
				return True
		</script>
</body>
</html>
'''

option_costs = {
		"option1": 3.42,
		"option2": 3.83,
		"option3": 4.40,
		"option4": 4.97,
		"option5": 5.54,
		"option6": 6.11,
		"option7": 6.52,
		"option8": 7.09,
		"option9": 7.66,
		"option10": 7.66,
		"option11": 8.23,
		"option12": 8.80,
		"option13": 10.35,
		"option14": 12.47,
		"option15": 13.04,
		"option16": 13.04,
		"option17": 10.92,
		"option18": 11.49,
		"option19": 14.18,
		"option20": 17.85,
		"option21": 17.85,
		"option22": 16.87,
		"option23": 19.56,
		"option24": 19.99
}

@app.route('/')
def index():
		return render_template_string(html_form)

@app.route('/create-checkout-session', methods=['POST'])
def create_checkout_session():
		selected_options = request.form.getlist('option')
		subscription_pack = request.form.get('subscription_pack')
		total_cost = sum(option_costs[option] for option in selected_options)

		if total_cost > 20:
				return f'Selection exceeds $20. Total cost: ${total_cost:.2f}. Please choose fewer options.'

		try:
				# Create Stripe Checkout session
				session = stripe.checkout.Session.create(
						payment_method_types=['card'],
						line_items=[
								{
										'price_data': {
												'currency': 'usd',
												'product_data': {
														'name': f'Subscription Pack {subscription_pack} with selected options',
												},
												'unit_amount': int(total_cost * 100),  # Convert dollars to cents
										},
										'quantity': 1,
								},
						],
						mode='payment',
						success_url=url_for('success', _external=True) + '?session_id={CHECKOUT_SESSION_ID}',
						cancel_url=url_for('index', _external=True),
				)
				return redirect(session.url, code=303)

		except Exception as e:
				return str(e)

@app.route('/success')
def success():
		return "Payment successful! Thank you for your purchase."

if __name__ == '__main__':
		app.run(debug=True, host='0.0.0.0', port=5000)
