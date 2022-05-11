import numpy as np
import matplotlib.pyplot as plt
import soundfile

sampleRate = 48000
seed = 9847265892659528763498

rng = np.random.default_rng(seed)
sig = 2.0 * rng.binomial(1, 0.5, sampleRate) - 1
soundfile.write(f"flip{seed}.wav", sig, sampleRate, subtype="FLOAT")
