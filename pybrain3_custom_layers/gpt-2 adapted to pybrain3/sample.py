# Dylan Kenneth Eliot & GPT-4-Plugins ( Alpha Edition )


"""
This produces the samples similar to the tensorflow equivalent. The only difference is it uses numpy to do the same calculation.



"""
import numpy as np


def top_k_logits(logits, k):
    if k == 0:
        return logits
    values = np.sort(logits, axis=-1)[:, -k:]
    min_values = values[:, -1, np.newaxis]
    return np.where(logits < min_values, np.full_like(logits, -1e10), logits)

def top_p_logits(logits, p):
    sorted_logits = np.sort(logits, axis=-1)[:, ::-1]
    cumulative_probs = np.cumsum(np.exp(sorted_logits) / np.sum(np.exp(sorted_logits), axis=-1, keepdims=True), axis=-1)
    indices_to_remove = cumulative_probs > p
    indices_to_remove[:, 1:] = indices_to_remove[:, :-1].copy()
    indices_to_remove[:, 0] = False
    sorted_indices = np.argsort(logits, axis=-1)[:, ::-1]
    rows = np.arange(len(logits))
    cols = sorted_indices[rows, np.argmax(indices_to_remove, axis=-1)]
    new_logits = logits.copy()
    new_logits[rows, cols] = -1e10
    return new_logits

def sample_sequence(hparams, length, start_token=None, batch_size=None, context=None, temperature=1, top_k=0, top_p=1):
    if start_token is None:
        assert context is not None, 'Specify exactly one of start_token and context!'
    else:
        assert context is None, 'Specify exactly one of start_token and context!'
        context = np.full((batch_size, 1), start_token)

    def step(hparams, tokens, past=None):
        # Assuming `model` is a global or passed reference to a model object that has been adapted for NumPy
        lm_output = model.model(hparams=hparams, X=tokens, past=past)

        logits = lm_output['logits'][:, :, :hparams.n_vocab]
        presents = lm_output['present']
        # No direct equivalent of set_shape in NumPy; shape management is implicit
        return {
            'logits': logits,
            'presents': presents,
        }

    def body(past, prev, output):
        next_outputs = step(hparams, prev, past=past)
        logits = next_outputs['logits'][:, -1, :] / temperature
        logits = top_k_logits(logits, k=top_k)  # Assume adapted for NumPy
        logits = top_p_logits(logits, p=top_p)  # Assume adapted for NumPy
        probs = np.exp(logits) / np.sum(np.exp(logits), axis=-1, keepdims=True)
        samples = np.array([np.random.choice(range(probs.shape[-1]), p=prob) for prob in probs], dtype=np.int32).reshape(batch_size, 1)
        return [
            np.concatenate([past, next_outputs['presents']], axis=-2) if past is not None else next_outputs['presents'],
            samples,
            np.concatenate([output, samples], axis=1)
        ]

    past, prev, output = body(None, context, context)

    for _ in range(length - 1):
        past, prev, output = body(past, prev, output)

    return output
