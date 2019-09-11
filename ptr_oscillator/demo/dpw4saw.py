import soundfile
from scipy import special

class DPW4Saw:
    def __init__(self, samplerate, frequency):
        N = 4

        self.tick = frequency / samplerate
        self.phase = 0
        self.y = [0] * N
        self.scale = 1 / (2 * self.tick)**(N - 1) / special.factorial(N)

        for _ in range(N):
            self.process()

    def process(self):
        x = 2 * self.phase - 1

        self.y[0] = self.y[1]
        self.y[1] = self.y[2]
        self.y[2] = self.y[3]
        self.y[3] = x * x * x * x - 2 * x * x

        output = self.y[0] - 3 * self.y[1] + 3 * self.y[2] - self.y[3]

        self.phase += self.tick
        if self.phase > 1:
            self.phase -= 1
        return self.scale * output


samplerate = 44100
frequency = 1000

osc = DPW4Saw(samplerate, frequency)
wav = [osc.process() for _ in range(samplerate)]

soundfile.write("snd/dpw4saw.wav", wav, samplerate, subtype="FLOAT")
