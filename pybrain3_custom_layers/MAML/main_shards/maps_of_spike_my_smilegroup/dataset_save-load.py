# Dylan Kenneth Eliot

from pybrain3.datasets.supervised import SupervisedDataSet

# 1. Create a simple supervised dataset
# Let's assume 2 inputs and 1 target
ds = SupervisedDataSet(2, 1)

# 2. Add some sample data
# Each addSample takes (input tuple, target tuple)
ds.addSample((0.0, 0.0), (0.0,))
ds.addSample((0.0, 1.0), (1.0,))
ds.addSample((1.0, 0.0), (1.0,))
ds.addSample((1.0, 1.0), (0.0,))

# 3. Save the dataset to a file
ds.saveToFile("my_dataset.ds")

print("Dataset saved as 'my_dataset.ds'")

# 4. (Optional) Load it back to verify
from pybrain3.datasets import SupervisedDataSet
loaded_ds = SupervisedDataSet.loadFromFile("my_dataset.ds")
print("Loaded dataset:")
for inpt, target in loaded_ds:
    print("Input:", inpt, "Target:", target)
