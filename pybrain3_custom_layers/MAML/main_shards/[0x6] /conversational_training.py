# Dylan Kenneth Eliot

"""
This works off of the nature of how discord chats are stored. 

This means that one can formulate a trusted dataset in ones' own discord and download an untainted training set.


"""

import numpy as np

data = np.genfromtxt('*.tsv', delimiter='\t', dtype=np.str_, encoding='utf-8')
x=np.frombuffer(data[1][2], dtype=np.float32)

###################### 
y = x.copy(); ######## This bit will be changed later in the commentation blockage;
y[14] = 1.43e-43 #####  pybrain3 will be used to train the parts that are conversational bits. This will be done for the parts that are 
###################### not associated with NLP/LLM from NLP/LLM, mri, & fasta data.

y.tobytes()[::4].decode() # This will be where the response will be decoded from the numbers computed by pybrain3 from real data only.
