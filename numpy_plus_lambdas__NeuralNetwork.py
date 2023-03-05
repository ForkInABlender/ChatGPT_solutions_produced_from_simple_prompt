import os
import numpy as np

# Load the Iris dataset
data = np.genfromtxt('iris.csv', delimiter=',', dtype='float')
X, y = data[:, :-1], data[:, -1]

# Preprocess the data
X = np.nan_to_num(X)  # Replace NaN values with 0
X = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0))  # Normalize the data to 0-1 range

# Split the data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Extract features
X_train_mean = np.mean(X_train, axis=0)
X_train_std = np.std(X_train, axis=0)
X_train_features = np.hstack((X_train_mean, X_train_std))

# Define the path to save/load the network weights
weights_path = "network_weights.npy"

if os.path.exists(weights_path):
    # Load the network weights if they exist
    theta = np.load(weights_path)
else:
    # Train a logistic regression model
    sigmoid = lambda z: 1 / (1 + np.exp(-z))

    theta = np.random.rand(X_train_features.shape[0])
    learning_rate = 0.01
    num_iterations = 1000

    # Use lambda function instead of for loop
    update_theta = lambda theta, _: theta - (learning_rate / y_train.size) * np.dot(X_train_features.T, (sigmoid(np.dot(X_train_features, theta)) - y_train))

    theta = np.fromiter(iter(lambda: update_theta(theta, None), theta), dtype=float, count=num_iterations)[-1]

    # Save the network weights
    np.save(weights_path, theta)

# Evaluate the model on the test set
X_test_mean = np.mean(X_test, axis=0)
X_test_std = np.std(X_test, axis=0)
X_test_features = np.hstack((X_test_mean, X_test_std))

z = np.dot(X_test_features, theta)
h = sigmoid(z)
y_pred = np.round(h)

accuracy = np.mean(y_pred == y_test)
print(f"Accuracy: {accuracy:.2f}")
