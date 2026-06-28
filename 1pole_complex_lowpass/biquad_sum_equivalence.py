import numpy as np


class ComplexLowpass:
    def __init__(self):
        self.x1 = 0
        self.y1 = 0

    def process(self, x0, cutoffNormalized, R):
        cut = np.exp(2j * np.pi * cutoffNormalized)
        res = np.exp(cut.imag * np.log(R))

        a = cut * res
        b = (1 - a) / 2

        self.y1 = b * (x0 + self.x1) + a * self.y1
        self.x1 = x0
        return self.y1


class BiquadFilter:
    def __init__(self):
        self.x1 = 0
        self.x2 = 0
        self.y1 = 0
        self.y2 = 0

    def process(self, x0, b0, b1, b2, a0, a1, a2):
        y0 = (b0 * x0 + b1 * self.x1 + b2 * self.x2 - a1 * self.y1 - a2 * self.y2) / a0
        self.x2 = self.x1
        self.x1 = x0
        self.y2 = self.y1
        self.y1 = y0
        return y0


np.random.seed(42)
inputs = np.random.randn(1000)

cutoffNormalized = 0.15
R = 0.9

clp = ComplexLowpass()
output_clp = np.array([clp.process(x, cutoffNormalized, R).real for x in inputs])

cut = np.exp(2j * np.pi * cutoffNormalized)
res = np.exp(cut.imag * np.log(R))
a = cut * res

alpha = (1 - abs(a) ** 2) / (1 + abs(a) ** 2)
cos_w0 = 2 * a.real / (1 + abs(a) ** 2)

a0 = 1 + alpha
a1 = -2 * cos_w0
a2 = 1 - alpha

b0_lpf = (1 - cos_w0) / 2
b1_lpf = 1 - cos_w0
b2_lpf = (1 - cos_w0) / 2

b0_bpf = alpha
b1_bpf = 0
b2_bpf = -alpha

lpf_filt = BiquadFilter()
bpf_filt = BiquadFilter()

output_biquads = []
for x in inputs:
    y_lpf = lpf_filt.process(x, b0_lpf, b1_lpf, b2_lpf, a0, a1, a2)
    y_bpf = bpf_filt.process(x, b0_bpf, b1_bpf, b2_bpf, a0, a1, a2)
    output_biquads.append(y_lpf + 0.5 * y_bpf)
output_biquads = np.array(output_biquads)

max_diff = np.max(np.abs(output_clp - output_biquads))
print(f"Max absolute difference: {max_diff:.2e}")
assert max_diff < 1e-12
