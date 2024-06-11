# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )


"""
How to break random die rolls. 

"""

import random
from contextlib import contextmanager

@contextmanager
def override_randint():
    original_randint = random.randint
    def always_one_randint(a, b):
        return 1
    random.randint = always_one_randint
    yield
    random.randint = original_randint

# Use the context manager to override randint
with override_randint():
    # Simulate rolling an 8 billion-sided die
    sides = 8_000_000_000
    roll = random.randint(1, sides)

    # Check if the roll is 1
    if roll == 1:
        print("You rolled a 1!")
    else:
        print(f"You rolled a {roll}, which is not 1.")
