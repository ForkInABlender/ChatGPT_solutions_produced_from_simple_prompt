# Dylan Kennth Eliot

"""
This is a dataset aimed at 3d input & 1d target datasets.

"""
from pybrain3.datasets.dataset import DataSet
import numpy as np

class Dim3DataSet(DataSet):
  """
  A custom dataset class to handle 3D input and 1D target data.
  """

  def __init__(self):
    """
    Initialize the dataset.
    """
    super().__init__()
    self.input_data = []
    self.target_data = []

  def add_sample(self, input_data, target_data):
    """
    Add a sample to the dataset.

    Args:
      input_data: A NumPy array representing the input data with shape (3,).
      target_data: A NumPy array representing the target data with shape (1,).
    """

    input_data = self._validate_data(input_data, (3,), "input")
    target_data = self._validate_data(target_data, (1,), "target")

    self.input_data.append(input_data)
    self.target_data.append(target_data)

  def _validate_data(self, data, expected_dim, name):
    """
    Validate the dimensions of the provided data.

    Args:
      data: A NumPy array to validate.
      expected_dim: A tuple of expected dimensions.
      name: A string name of the data being validated ('input' or 'target').

    Returns:
      The validated data.

    Raises:
      ValueError: If the data is not a NumPy array or has incorrect dimensions.
    """

    if not isinstance(data, np.ndarray):
      raise ValueError(f"{name.capitalize()} data must be a numpy array.")

    if data.shape != expected_dim:
      raise ValueError(f"{name.capitalize()} data shape {data.shape} does not match expected dimensions {expected_dim}.")

    return data

  def get_samples(self):
    """
    Retrieve all samples in the dataset.

    Returns:
      A list of tuples containing (input_data, target_data).
    """

    return list(zip(self.input_data, self.target_data))

  def __len__(self):
    """
    Get the number of samples in the dataset.

    Returns:
      An integer count of the samples.
    """

    return len(self.input_data)

  def __getitem__(self, index):
    """
    Get a specific sample by index.

    Args:
      index: The integer index of the sample.

    Returns:
      A tuple containing (input_data, target_data).
    """

    return self.input_data[index], self.target_data[index]
