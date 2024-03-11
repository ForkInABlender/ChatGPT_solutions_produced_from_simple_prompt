# Dylan Kenneth Eliot & GPT-4-Plugins ( Alpha Edition )

"""

This now loads the numpy version of gpt-2 as before but from unconditional to conditional interactive like state. But without the tensorflow overhead. 

"""


import fire
import json
import os
import numpy as np

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

def interact_model(
    model_name='124M',
    seed=None,
    nsamples=1,
    batch_size=1,
    length=None,
    temperature=1,
    top_k=0,
    top_p=1,
    models_dir='models',
):
    """
    Interactively run the model using NumPy/SciPy logic
    """
    models_dir = os.path.expanduser(os.path.expandvars(models_dir))
    assert nsamples % batch_size == 0

    enc = encoder.get_encoder(model_name, models_dir)
    hparams = model.default_hparams()
    with open(os.path.join(models_dir, model_name, 'hparams.json')) as f:
        hparams.override_from_dict(json.load(f))

    if length is None:
        length = hparams.n_ctx // 2
    elif length > hparams.n_ctx:
        raise ValueError("Can't get samples longer than window size: %s" % hparams.n_ctx)

    np.random.seed(seed)

    # Load your model weights here into a Python object
    my_model_weights = load_my_model(model_name, models_dir)

    while True:
        raw_text = input("Model prompt >>> ")
        while not raw_text:
            print('Prompt should not be empty!')
            raw_text = input("Model prompt >>> ")
        context_tokens = enc.encode(raw_text)
        generated = 0
        for _ in range(nsamples // batch_size):
            # Assuming `sample_sequence` is adapted to use NumPy for generating a sample sequence
            out = sample.sample_sequence(
                model_weights=my_model_weights,  # Pass the loaded model weights
                hparams=hparams,
                length=length + len(context_tokens),
                context=context_tokens,
                batch_size=batch_size,
                temperature=temperature,
                top_k=top_k,
                top_p=top_p
            )

            for i in range(batch_size):
                generated += 1
                text = enc.decode(out[i][len(context_tokens):])  # Skip the context part of the output
                print("=" * 40 + " SAMPLE " + str(generated) + " " + "=" * 40)
                print(text)
            print("=" * 80)
if __name__ == '__main__':
    fire.Fire(interact_model)
