# Dylan Kenneth Eliot & GPT-4-Plugins ( Alpha Edition )

"""
This uses numpy & scipy to do what is normally done with tensorflow.

There is also a way to save at the bottom of the file as well.

"""


import fire
import json
import os
import numpy as np
import scipy.stats as stats

# Assuming `model`, `sample`, and `encoder` modules are adapted for NumPy/SciPy usage
import model, sample, encoder

def load_my_model(model_name, models_dir):
    """
    Load model weights from a .npz file into a dictionary of NumPy arrays.
    """
    weights_path = os.path.join(models_dir, model_name, 'model_weights_dict.npz')
    with np.load(weights_path) as data:
        model_weights = {key: data[key] for key in data}
    return model_weights

def sample_model(
    model_name='124M',
    seed=None,
    nsamples=0,
    batch_size=1,
    length=None,
    temperature=1,
    top_k=0,
    top_p=1,
    models_dir='models',
):
    """
    Run the sample_model with NumPy/SciPy logic
    """
    models_dir = os.path.expanduser(os.path.expandvars(models_dir))
    enc = encoder.get_encoder(model_name, models_dir)
    hparams = model.default_hparams()
    with open(os.path.join(models_dir, model_name, 'hparams.json')) as f:
        hparams.override_from_dict(json.load(f))

    if length is None:
        length = hparams.n_ctx
    elif length > hparams.n_ctx:
        raise ValueError("Can't get samples longer than window size: %s" % hparams.n_ctx)

    np.random.seed(seed)

    # Load your model weights here into a Python object
    my_model_weights = load_my_model(model_name, models_dir)

    generated = 0
    while nsamples == 0 or generated < nsamples:
        # Assuming `sample_sequence` is adapted to use NumPy for generating a sample sequence
        out = sample.sample_sequence(
            model_weights=my_model_weights,  # Pass the loaded model weights
            hparams=hparams,
            length=length,
            start_token=enc.encoder[''],
            batch_size=batch_size,
            temperature=temperature,
            top_k=top_k,
            top_p=top_p
        )

        for i in range(batch_size):
            generated += batch_size
            text = enc.decode(out[i])
            print("=" * 40 + " SAMPLE " + str(generated) + " " + "=" * 40)
            print(text)
if __name__ == '__main__':
    fire.Fire(sample_model)

"""
def save_my_model(model_weights, model_name, models_dir):
    # Ensure the directory exists
    model_dir = os.path.join(models_dir, model_name)
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)
    
    # Construct the path to the weights file
    weights_path = os.path.join(model_dir, 'model_weights_dict.npz')
    
    # Save the weights
    np.savez_compressed(weights_path, **model_weights)



"""
