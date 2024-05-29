from flask import Flask, render_template, request, redirect, url_for, session
import stripe
import os
import sqlite3
from datetime import datetime, timedelta

app = Flask(__name__)
app.secret_key = os.urandom(24)  # Secret key for session management

# Stripe configuration
stripe_keys = {
    'secret_key': 'your_stripe_secret_key',
    'publishable_key': 'your_stripe_publishable_key'
}

stripe.api_key = stripe_keys['secret_key']

# Sample product data
products = [
    {'id': 1, 'name': 'Product 1', 'price': 1000},  # Price in cents
    {'id': 2, 'name': 'Product 2', 'price': 2000},
    {'id': 3, 'name': 'Product 3', 'price': 3000}
]

# Check if the user needs to pay based on the 30-day timeout
def needs_payment(email):
    conn = sqlite3.connect('payments.db')
    c = conn.cursor()
    c.execute("SELECT timestamp FROM payments WHERE email = ? ORDER BY timestamp DESC LIMIT 1", (email,))
    row = c.fetchone()
    conn.close()
    if row:
        last_payment = datetime.strptime(row[0], '%Y-%m-%d %H:%M:%S')
        if datetime.now() - last_payment < timedelta(days=30):
            return False
    return True

# Record a payment in the database
def record_payment(email):
    conn = sqlite3.connect('payments.db')
    c = conn.cursor()
    c.execute("INSERT INTO payments (email, timestamp) VALUES (?, ?)", (email, datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    conn.commit()
    conn.close()

@app.route('/')
def index():
    return render_template('index.html', products=products)

@app.route('/add_to_cart/<int:product_id>')
def add_to_cart(product_id):
    product = next((item for item in products if item['id'] == product_id), None)
    if product:
        cart = session.get('cart', [])
        cart.append(product)
        session['cart'] = cart
    return redirect(url_for('index'))

@app.route('/cart')
def cart():
    cart = session.get('cart', [])
    total_amount = sum(item['price'] for item in cart)
    return render_template('cart.html', cart=cart, total_amount=total_amount, key=stripe_keys['publishable_key'])

@app.route('/charge', methods=['POST'])
def charge():
    email = request.form['stripeEmail']
    if needs_payment(email):
        cart = session.get('cart', [])
        total_amount = sum(item['price'] for item in cart)

        customer = stripe.Customer.create(
            email=email,
            source=request.form['stripeToken']
        )

        charge = stripe.Charge.create(
            customer=customer.id,
            amount=total_amount,
            currency='usd',
            description='POS Charge'
        )

        record_payment(email)
        session.pop('cart', None)  # Clear the cart
        return redirect(url_for('success', amount=total_amount))
    else:
        return redirect(url_for('index'))

@app.route('/success/<int:amount>')
def success(amount):
    return render_template('success.html', amount=amount)

if __name__ == '__main__':
    app.run(debug=True)
