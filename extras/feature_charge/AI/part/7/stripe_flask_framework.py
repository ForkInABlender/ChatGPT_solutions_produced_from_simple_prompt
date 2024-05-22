# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
Welcome to the basics of what openai looks like with the final part.

This will be used for features available later on as one decides on more formidable layout of charge. Each feature will be billed separate but the 
 base bill will be based on what they use.

"""


from flask import Flask, render_template, request, redirect, url_for, session, flash
from flask_sqlalchemy import SQLAlchemy
from flask_login import LoginManager, UserMixin, login_user, login_required, logout_user, current_user
import stripe
import os
from datetime import datetime, timedelta

app = Flask(__name__)
app.secret_key = os.urandom(24)

# Database setup
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///payments.db'
db = SQLAlchemy(app)

# Flask-Login setup
login_manager = LoginManager()
login_manager.init_app(app)
login_manager.login_view = 'login'

# Stripe configuration
stripe_keys = {
    'secret_key': 'your_stripe_secret_key',
    'publishable_key': 'your_stripe_publishable_key'
}

stripe.api_key = stripe_keys['secret_key']

# User model
class User(UserMixin, db.Model):
    id = db.Column(db.Integer, primary_key=True)
    email = db.Column(db.String(150), unique=True, nullable=False)
    password = db.Column(db.String(150), nullable=False)
    last_payment = db.Column(db.DateTime, nullable=True)
    features = db.relationship('FeatureUsage', backref='user', lazy=True)

# Payment model
class Payment(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.Integer, db.ForeignKey('user.id'), nullable=False)
    timestamp = db.Column(db.DateTime, nullable=False)
    amount = db.Column(db.Integer, nullable=False)

# Feature model
class Feature(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(100), nullable=False)
    cost = db.Column(db.Integer, nullable=False)

# Feature usage model
class FeatureUsage(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.Integer, db.ForeignKey('user.id'), nullable=False)
    feature_id = db.Column(db.Integer, db.ForeignKey('feature.id'), nullable=False)
    timestamp = db.Column(db.DateTime, nullable=False, default=datetime.utcnow)
    billed = db.Column(db.Boolean, default=False)

# Create tables
with app.app_context():
    db.create_all()

@login_manager.user_loader
def load_user(user_id):
    return User.query.get(int(user_id))

@app.route('/login', methods=['GET', 'POST'])
def login():
    if request.method == 'POST':
        email = request.form['email']
        password = request.form['password']
        user = User.query.filter_by(email=email).first()
        if user and user.password == password:
            login_user(user)
            return redirect(url_for('index'))
        else:
            flash('Login failed. Check your email and password.')
    return render_template('login.html')

@app.route('/logout')
@login_required
def logout():
    logout_user()
    return redirect(url_for('login'))

@app.route('/')
@login_required
def index():
    return render_template('index.html')

@app.route('/use_feature/<int:feature_id>')
@login_required
def use_feature(feature_id):
    feature = Feature.query.get(feature_id)
    if feature:
        usage = FeatureUsage(user_id=current_user.id, feature_id=feature_id)
        db.session.add(usage)
        db.session.commit()
        flash(f'You used {feature.name}. Cost: ${feature.cost / 100}.')
        return redirect(url_for('index'))
    else:
        flash('Feature not found.')
        return redirect(url_for('index'))

@app.route('/charge', methods=['POST'])
@login_required
def charge():
    email = current_user.email
    usages = FeatureUsage.query.filter_by(user_id=current_user.id, billed=False).all()
    total_amount = sum(usage.feature.cost for usage in usages)

    if total_amount > 0:
        customer = stripe.Customer.create(
            email=email,
            source=request.form['stripeToken']
        )

        charge = stripe.Charge.create(
            customer=customer.id,
            amount=total_amount,
            currency='usd',
            description='Feature Usage Charge'
        )

        payment = Payment(user_id=current_user.id, timestamp=datetime.utcnow(), amount=total_amount)
        db.session.add(payment)
        for usage in usages:
            usage.billed = True
        db.session.commit()

        flash(f'You have been charged ${total_amount / 100}.')
        return redirect(url_for('index'))
    else:
        flash('No charges to process.')
        return redirect(url_for('index'))

@app.route('/success/<int:amount>')
@login_required
def success(amount):
    return render_template('success.html', amount=amount)

if __name__ == '__main__':
    app.run(debug=True)
