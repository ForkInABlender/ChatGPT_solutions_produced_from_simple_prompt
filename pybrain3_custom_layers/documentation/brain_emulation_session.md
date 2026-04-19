# Brain Emulation Pipeline — Session Notes
**Date:** 2026-04-18

---

## Stack

- **NumPy + SciPy** — numerical foundation, signal processing, spike train analysis
- **RDKit** — molecular descriptors (molwt, logP) from SMILES; drives neuron resting potential
- **PyBrain3** — neural network layers, backprop trainer, custom neuron layers
- **pandas + odfpy** — ODS spreadsheet storage/retrieval per brain region

---

## Core Neuron Model

Reference: [`AI_neurons.py`](https://raw.githubusercontent.com/ForkInABlender/ChatGPT_solutions_produced_from_simple_prompt/refs/heads/2026_03/pybrain3_custom_layers/AI_neurons.py)

- RK4 `Integrator` class drives membrane voltage: `dv/dt = (molwt + logP - v) / tau`
- Molecular weight + logP from SMILES directly set neuron resting potential (neurochemical coupling)
- Spike detection: `v >= v_th` → record spike, reset to `v_reset`
- PyBrain3 network trained on spike states → predicts molecular weight as supervised signal

---

## Brain Region Model

Reference: [`Brain_emulator.py`](https://raw.githubusercontent.com/ForkInABlender/ChatGPT_solutions_produced_from_simple_prompt/refs/heads/2026_03/pybrain3_custom_layers/Brain_emulator.py)  
Reference: [`main__GPT_Model.py`](https://raw.githubusercontent.com/ForkInABlender/ChatGPT_solutions_produced_from_simple_prompt/refs/heads/2026_03/pybrain3_custom_layers/main__GPT_Model.py)

### Regions (17 total)

| # | Region | Dims | RAM | ODS on disk |
|---|---|---|---|---|
| 1 | Anterior Cingulate Cortex | 100→250→50 recurrent | ~0.3 MB | ~0.1 MB |
| 2 | Mirror Neuron System | 500→250→50 recurrent | ~0.6 MB | ~0.2 MB |
| 3 | Insula | 500→1→50 | ~0.1 MB | ~0.04 MB |
| 4 | Frontal Lobe | 100→1→50257 | ~20 MB | ~1.5 MB |
| 5 | Prefrontal Cortex | 230→1→100 | ~0.1 MB | ~0.04 MB |
| 6 | Cerebral Cortex | 151386→1→50257 | ~61 MB | ~5 MB |
| 7 | Temporal Lobe | 50→1→250 | ~0.05 MB | ~0.02 MB |
| 8 | Visual Cortex | 120→1→80 | ~0.04 MB | ~0.02 MB |
| 9 | Broca's Area | 150→1→75 | ~0.05 MB | ~0.02 MB |
| 10 | Wernicke's Area | 75→1→50257 | ~15 MB | ~1.2 MB |
| 11 | Cerebellum | 75→1→70 | ~0.02 MB | ~0.01 MB |
| 12 | Brainstem | 120→1→60 | ~0.03 MB | ~0.01 MB |
| 13 | Hippocampus | 250→1→50 LSTM | ~0.1 MB | ~0.04 MB |
| 14 | Thalamus | 400→1→120 | ~0.2 MB | ~0.08 MB |
| 15 | Amygdala | 120→1→10 | ~0.005 MB | ~0.002 MB |
| 16 | Hypothalamus | 200→1→20 | ~0.02 MB | ~0.01 MB |
| 17 | GPT (LLM sandwich) | 50257→128×96 blocks | ~410 MB | ~40 MB |
| | **TOTAL** | | **~508 MB** | **~48 MB** |

### Notable absences from full brain atlas
Parietal lobe, occipital lobe, basal ganglia, corpus callosum, olfactory bulb, motor cortex, somatosensory cortex.

---

## Pipeline Architecture

```
.nii MRI
    ↓ nibabel → 3D activation volume → threshold → activation mask
    ↓ directed graph (nodes=regions, edges=activation-ordered connections)
         ↓
    Per-node: SMILES of dominant neurotransmitter/receptor
         ↓ RDKit → molwt, logP
         ↓ Integrator(y0, deriv, t, dt) → spike pattern
         ↓
    Wernicke subgraph → spike encoding → latent vector
         ↓
    [GPT language model — sandwiched]
         ↓
    latent vector → spike decoding
         ↓
    Broca subgraph → motor/speech output spikes
```

---

## Neurochemical Exchange Data Per Region

| # | Region | Primary neurotransmitters | Key receptors (FASTA) | MRI size | FASTA size | ODS on disk |
|---|---|---|---|---|---|---|
| 1 | ACC | Glutamate, GABA, Dopamine | NMDA, GABA-A, D2 | ~5 MB | ~0.3 MB | ~0.8 MB |
| 2 | Mirror Neuron | Glutamate, Dopamine | NMDA, D1 | ~8 MB | ~0.2 MB | ~1.2 MB |
| 3 | Insula | Serotonin, Acetylcholine | 5-HT2A, nAChR | ~6 MB | ~0.2 MB | ~0.9 MB |
| 4 | Frontal Lobe | Dopamine, Glutamate | D1, NMDA, AMPA | ~15 MB | ~0.4 MB | ~2.2 MB |
| 5 | Prefrontal Cortex | Dopamine, Norepinephrine | D1, α2A | ~12 MB | ~0.3 MB | ~1.8 MB |
| 6 | Cerebral Cortex | Glutamate, GABA | AMPA, GABA-A | ~50 MB | ~1.2 MB | ~7.5 MB |
| 7 | Temporal Lobe | Glutamate, Acetylcholine | NMDA, mAChR | ~10 MB | ~0.3 MB | ~1.5 MB |
| 8 | Visual Cortex | Glutamate, GABA | AMPA, GABA-B | ~30 MB | ~0.5 MB | ~4.5 MB |
| 9 | Broca's Area | Glutamate, GABA | AMPA, GABA-A | ~5 MB | ~0.2 MB | ~0.8 MB |
| 10 | Wernicke's Area | Glutamate, GABA | NMDA, AMPA | ~5 MB | ~0.2 MB | ~0.8 MB |
| 11 | Cerebellum | GABA, Glutamate | GABA-A, mGluR | ~20 MB | ~0.4 MB | ~3.0 MB |
| 12 | Brainstem | Serotonin, Norepinephrine | 5-HT1A, α1 | ~8 MB | ~0.2 MB | ~1.2 MB |
| 13 | Hippocampus | Glutamate, Acetylcholine | NMDA, mAChR | ~10 MB | ~0.3 MB | ~1.5 MB |
| 14 | Thalamus | Glutamate, GABA | AMPA, GABA-B | ~12 MB | ~0.3 MB | ~1.8 MB |
| 15 | Amygdala | GABA, Dopamine, Serotonin | GABA-A, D2, 5-HT2A | ~6 MB | ~0.2 MB | ~0.9 MB |
| 16 | Hypothalamus | Dopamine, Oxytocin, AVP | D2, OXTR, V1aR | ~5 MB | ~0.3 MB | ~0.8 MB |
| 17 | GPT (LLM) | n/a | n/a | n/a | n/a | n/a |
| | **TOTAL** | | | **~207 MB** | **~5.5 MB** | **~31 MB** |

---

## Storage Strategy

### Format
- **Serialization:** JSON (human-inspectable, not pickle)
- **Encoding:** base64
- **Cipher:** gridshift
- **Container:** ODS via pandas + odfpy
- **Cell limit:** 50,000 chars/cell, 1,000×1,000 cells/sheet

### Save/Load pattern
```python
import pandas as pd
import base64, json

def save_region(region_id, state_dict):
    encoded = base64.b64encode(json.dumps(state_dict).encode()).decode()
    chunks = [encoded[i:i+50000] for i in range(0, len(encoded), 50000)]
    pd.DataFrame({"data": chunks}).to_excel(f"region_{region_id}.ods", engine="odf", index=False)

def load_region(region_id):
    df = pd.read_excel(f"region_{region_id}.ods", engine="odf")
    encoded = "".join(df["data"].tolist())
    return json.loads(base64.b64decode(encoded).decode())
```

### One ODS file per region — lazy load on activation, evict on silence

---

## Total Memory Summary

| Component | RAM | ODS on disk | Compacted ODS |
|---|---|---|---|
| Network weights (17 regions) | ~508 MB | ~48 MB | ~38 MB |
| MRI activation data | ~207 MB | ~25 MB | ~12 MB |
| FASTA + SMILES | ~5.6 MB | ~1 MB | ~1.5 MB |
| **Total** | **~721 MB** | **~74 MB** | **~52 MB** |

### With lazy loading
- gpt_net always loaded: ~410 MB
- 2–3 active regions: ~50 MB
- Neurochemical data: ~20 MB
- **Working RAM: ~480 MB**

### Compression pipeline
Raw → positional encoder (occurrence dict, ~20–80% reduction by data type) → JSON → base64 (+33%) → gridshift → ODS ZIP (~65% reduction on XML)
Net: **~28% smaller than uncompacted ODS on disk**

---

## Pyodide / In-Browser Compatibility

| Component | Status |
|---|---|
| NumPy, SciPy | ✓ |
| PyBrain3 | ✓ |
| JSON, base64, zipfile | ✓ |
| odfpy | ✓ (pure Python) |
| Spike/Integrator class | ✓ |
| RDKit | swap for `rdkit.js` |
| nibabel | not in Pyodide — pre-process .nii data externally |
| Filesystem | use Pyodide virtual FS or IndexedDB |

**Working RAM (~480 MB) fits within browser tab limit (2–4 GB). Runs in-browser via Pyodide with lazy loading.**

---

## Memory Reduction Journey

| Step | Working RAM |
|---|---|
| Full neuron-for-neuron fidelity (86B neurons) | ~500 TB |
| Region-population abstraction (1M neurons) | ~15–20 GB |
| Lazy cell loading (50–100 active regions) | ~1–2 GB |
| Sharded ODS per region + actual model scale | ~480 MB |
| **Total reduction** | **~99.9999%** |
