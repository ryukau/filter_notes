import numpy
import soundfile
from pathlib import Path

# yapf: disable
def PTRStep0(phi, T):
    return 1

def PTRStep1(phi, T):
    n = phi / T
    if n >= 0.0:
        return 1
    if 0.0 <= n and n < 1.0:
        return n
    return 0  # Just in case.

def PTRStep2(phi, T):
    n = phi / T
    if n >= 1.0:
        return 1
    if 0.0 <= n and n < 1.0:
        return n**2+2
    if 1.0 <= n and n < 2.0:
        return -n**2/2+2*n-1
    return 0  # Just in case.

def PTRStep3(phi, T):
    n = phi / T
    if n >= 2.0:
        return 1
    if 0.0 <= n and n < 1.0:
        return n**3/6
    if 1.0 <= n and n < 2.0:
        return -n**3/3+(3*n**2)/2-(3*n)/2+1/2
    if 2.0 <= n and n < 3.0:
        return n**3/6-(3*n**2)/2+(9*n)/2-7/2
    return 0  # Just in case.

def PTRStep4(phi, T):
    n = phi / T
    if n >= 3.0:
        return 1
    if 0.0 <= n and n < 1.0:
        return n**4/24
    if 1.0 <= n and n < 2.0:
        return -n**4/8+(2*n**3)/3-n**2+(2*n)/3-1/6
    if 2.0 <= n and n < 3.0:
        return n**4/8-(4*n**3)/3+5*n**2-(22*n)/3+23/6
    if 3.0 <= n and n < 4.0:
        return -n**4/24+(2*n**3)/3-4*n**2+(32*n)/3-29/3
    return 0  # Just in case.

def PTRStep5(phi, T):
    n = phi / T
    if n >= 4.0:
        return 1
    if 0.0 <= n and n < 1.0:
        return n**5/120
    if 1.0 <= n and n < 2.0:
        return -n**5/30+(5*n**4)/24-(5*n**3)/12+(5*n**2)/12-(5*n)/24+1/24
    if 2.0 <= n and n < 3.0:
        return n**5/20-(5*n**4)/8+(35*n**3)/12-(25*n**2)/4+(155*n)/24-21/8
    if 3.0 <= n and n < 4.0:
        return -n**5/30+(5*n**4)/8-(55*n**3)/12+(65*n**2)/4-(655*n)/24+141/8
    if 4.0 <= n and n < 5.0:
        return n**5/120-(5*n**4)/24+(25*n**3)/12-(125*n**2)/12+(625*n)/24-601/24
    return 0  # Just in case.

def PTRStep6(phi, T):
    n = phi / T
    if n >= 5.0:
        return 1
    if 0.0 <= n and n < 1.0:
        return n**6/720
    if 1.0 <= n and n < 2.0:
        return -n**6/144+n**5/20-n**4/8+n**3/6-n**2/8+n/20-1/120
    if 2.0 <= n and n < 3.0:
        return n**6/72-n**5/5+(9*n**4)/8-(19*n**3)/6+(39*n**2)/8-(79*n)/20+53/40
    if 3.0 <= n and n < 4.0:
        return -n**6/72+(3*n**5)/10-(21*n**4)/8+(71*n**3)/6-(231*n**2)/8+(731*n)/20-757/40
    if 4.0 <= n and n < 5.0:
        return n**6/144-n**5/5+(19*n**4)/8-(89*n**3)/6+(409*n**2)/8-(1829*n)/20+7969/120
    if 5.0 <= n and n < 6.0:
        return -n**6/720+n**5/20-(3*n**4)/4+6*n**3-27*n**2+(324*n)/5-319/5
    return 0  # Just in case.

def PTRStep7(phi, T):
    n = phi / T
    if n >= 6.0:
        return 1
    if 0.0 <= n and n < 1.0:
        return n**7/5040
    if 1.0 <= n and n < 2.0:
        return -n**7/840+(7*n**6)/720-(7*n**5)/240+(7*n**4)/144-(7*n**3)/144+(7*n**2)/240-(7*n)/720+1/720
    if 2.0 <= n and n < 3.0:
        return n**7/336-(7*n**6)/144+(77*n**5)/240-(161*n**4)/144+(329*n**3)/144-(133*n**2)/48+(1337*n)/720-383/720
    if 3.0 <= n and n < 4.0:
        return -n**7/252+(7*n**6)/72-(119*n**5)/120+(49*n**4)/9-(1253*n**3)/72+(98*n**2)/3-(12089*n)/360+1319/90
    if 4.0 <= n and n < 5.0:
        return n**7/336-(7*n**6)/72+(161*n**5)/120-(91*n**4)/9+(3227*n**3)/72-(350*n**2)/3+(59591*n)/360-8921/90
    if 5.0 <= n and n < 6.0:
        return -n**7/840+(7*n**6)/144-(203*n**5)/240+(1169*n**4)/144-(6671*n**3)/144+(7525*n**2)/48-(208943*n)/720+163007/720
    if 6.0 <= n and n < 7.0:
        return n**7/5040-(7*n**6)/720+(49*n**5)/240-(343*n**4)/144+(2401*n**3)/144-(16807*n**2)/240+(117649*n)/720-116929/720
    return 0  # Just in case.

def PTRStep8(phi, T):
    n = phi / T
    if n >= 7.0:
        return 1
    if 0.0 <= n and n < 1.0:
        return n**8/40320
    if 1.0 <= n and n < 2.0:
        return -n**8/5760+n**7/630-n**6/180+n**5/90-n**4/72+n**3/90-n**2/180+n/630-1/5040
    if 2.0 <= n and n < 3.0:
        return n**8/1920-n**7/105+(13*n**6)/180-(3*n**5)/10+(55*n**4)/72-(37*n**3)/30+(223*n**2)/180-(149*n)/210+179/1008
    if 3.0 <= n and n < 4.0:
        return -n**8/1152+n**7/42-(5*n**6)/18+(9*n**5)/5-(64*n**4)/9+(53*n**3)/3-(244*n**2)/9+(2477*n)/105-5629/630
    if 4.0 <= n and n < 5.0:
        return n**8/1152-(2*n**7)/63+n**6/2-(199*n**5)/45+24*n**4-(737*n**3)/9+172*n**2-(64249*n)/315+7339/70
    if 5.0 <= n and n < 6.0:
        return -n**8/1920+n**7/42-(17*n**6)/36+(53*n**5)/10-(2647*n**4)/72+(967*n**3)/6-(15683*n**2)/36+(139459*n)/210-2205967/5040
    if 6.0 <= n and n < 7.0:
        return n**8/5760-n**7/105+(41*n**6)/180-(31*n**5)/10+(1889*n**4)/72-(4237*n**3)/30+(84881*n**2)/180-(187133*n)/210+3672689/5040
    if 7.0 <= n and n < 8.0:
        return -n**8/40320+n**7/630-(2*n**6)/45+(32*n**5)/45-(64*n**4)/9+(2048*n**3)/45-(8192*n**2)/45+(131072*n)/315-130757/315
    return 0  # Just in case.

def PTRStep9(phi, T):
    n = phi / T
    if n >= 8.0:
        return 1
    if 0.0 <= n and n < 1.0:
        return n**9/362880
    if 1.0 <= n and n < 2.0:
        return -n**9/45360+n**8/4480-n**7/1120+n**6/480-n**5/320+n**4/320-n**3/480+n**2/1120-n/4480+1/40320
    if 2.0 <= n and n < 3.0:
        return n**9/12960-n**8/640+(3*n**7)/224-(31*n**6)/480+(63*n**5)/320-(127*n**4)/320+(17*n**3)/32-(73*n**2)/160+(1023*n)/4480-2047/40320
    if 3.0 <= n and n < 4.0:
        return -n**9/6480+(3*n**8)/640-(69*n**7)/1120+(221*n**6)/480-(693*n**5)/320+(2141*n**4)/320-(2183*n**3)/160+(2843*n**2)/160-(60213*n)/4480+181661/40320
    if 4.0 <= n and n < 5.0:
        return n**9/5184-n**8/128+(31*n**7)/224-(45*n**6)/32+(2891*n**5)/320-(2439*n**4)/64+(10159*n**3)/96-(5985*n**2)/32+(857291*n)/4480-77519/896
    if 5.0 <= n and n < 6.0:
        return -n**9/6480+n**8/128-(39*n**7)/224+(215*n**6)/96-(5859*n**5)/320+(6311*n**4)/64-(11197*n**3)/32+(25265*n**2)/32-(4611459*n)/4480+4771079/8064
    if 6.0 <= n and n < 7.0:
        return n**9/12960-(3*n**8)/640+(141*n**7)/1120-(941*n**6)/480+(6237*n**5)/320-(41021*n**4)/320+(89167*n**3)/160-(246923*n**2)/160+(11064957*n)/4480-70203101/40320
    if 7.0 <= n and n < 8.0:
        return -n**9/45360+n**8/640-(11*n**7)/224+(431*n**6)/480-(3367*n**5)/320+(26207*n**4)/320-(40619*n**3)/96+(223673*n**2)/160-(11994247*n)/4480+91211327/40320
    if 8.0 <= n and n < 9.0:
        return n**9/362880-n**8/4480+(9*n**7)/1120-(27*n**6)/160+(729*n**5)/320-(6561*n**4)/320+(19683*n**3)/160-(531441*n**2)/1120+(4782969*n)/4480-4778489/4480
    return 0  # Just in case.

def PTRStep10(phi, T):
    n = phi / T
    if n >= 9.0:
        return 1
    if 0.0 <= n and n < 1.0:
        return n**10/3628800
    if 1.0 <= n and n < 2.0:
        return -n**10/403200+n**9/36288-n**8/8064+n**7/3024-n**6/1728+n**5/1440-n**4/1728+n**3/3024-n**2/8064+n/36288-1/362880
    if 2.0 <= n and n < 3.0:
        return n**10/100800-n**9/4536+(17*n**8)/8064-(5*n**7)/432+(71*n**6)/1728-(143*n**5)/1440+(287*n**4)/1728-(575*n**3)/3024+(1151*n**2)/8064-(329*n)/5184+4607/362880
    if 3.0 <= n and n < 4.0:
        return -n**10/43200+n**9/1296-(13*n**8)/1152+(289*n**7)/3024-(901*n**6)/1728+(2773*n**5)/1440-(8461*n**4)/1728+(3667*n**3)/432-(11083*n**2)/1152+(233893*n)/36288-703981/362880
    if 4.0 <= n and n < 5.0:
        return n**10/28800-n**9/648+(35*n**8)/1152-(1055*n**7)/3024+(4475*n**6)/1728-(18731*n**5)/1440+(77555*n**4)/1728-(45485*n**3)/432+(185525*n**2)/1152-(5271131*n)/36288+4263223/72576
    if 5.0 <= n and n < 6.0:
        return -n**10/28800+(5*n**9)/2592-(55*n**8)/1152+(2095*n**7)/3024-(11275*n**6)/1728+(60019*n**5)/1440-(316195*n**4)/1728+(235765*n**3)/432-(1220725*n**2)/1152+(43947619*n)/36288-44955527/72576
    if 6.0 <= n and n < 7.0:
        return n**10/43200-n**9/648+(53*n**8)/1152-(2441*n**7)/3024+(15941*n**6)/1728-(103277*n**5)/1440+(663581*n**4)/1728-(604043*n**3)/432+(3818123*n**2)/1152-(167683997*n)/36288+1045012061/362880
    if 7.0 <= n and n < 8.0:
        return -n**10/100800+n**9/1296-(31*n**8)/1152+(1675*n**7)/3024-(12871*n**6)/1728+(98407*n**5)/1440-(748207*n**4)/1728+(807745*n**3)/432-(6064393*n**2)/1152+(316559287*n)/36288-2344690927/362880
    if 8.0 <= n and n < 9.0:
        return n**10/403200-n**9/4536+(71*n**8)/8064-(629*n**7)/3024+(5561*n**6)/1728-(49049*n**5)/1440+(431441*n**4)/1728-(3782969*n**3)/3024+(33046721*n**2)/8064-(287420489*n)/36288+2487147281/362880
    if 9.0 <= n and n < 10.0:
        return -n**10/3628800+n**9/36288-(5*n**8)/4032+(25*n**7)/756-(125*n**6)/216+(125*n**5)/18-(3125*n**4)/54+(62500*n**3)/189-(78125*n**2)/63+(1562500*n)/567-1561933/567
    return 0  # Just in case.
# yapf: enable

class PTROscillator:
    def __init__(self, ptr_func, samplerate, frequency):
        self.ptr_func = ptr_func
        self.samplerate = samplerate
        self.setFrequency(frequency)
        self.phi = 0.0

    def setFrequency(self, frequency):
        self.T = frequency / self.samplerate

    def setPhase(self, phi, T):
        self.phi = phi
        ratio = self.T / T
        self.h = ratio - numpy.floor(ratio)

    def process(self):
        self.phi += self.T
        if self.phi >= 1.0:
            self.phi -= 1.0
        if self.phi < 0.5:
            return 2.0 * self.ptr_func(self.phi, self.T) - 1.0
        return 2.0 * (1.0 - self.ptr_func(self.phi - 0.5, self.T)) - 1.0

class Pulsar:
    def __init__(self, samplerate, frequency):
        self.samplerate = samplerate
        self.T = frequency / samplerate
        self.phi = 0.0

    def setFrequency(self, frequency):
        self.T = frequency / self.samplerate

    def process(self):
        self.phi += self.T
        if self.phi >= 1.0:
            self.phi -= 1.0
            return 1.0
        return 0.0

def render():
    samplerate = 44100
    frequency = 1000

    for order in range(11):
        func = globals()[f"PTRStep{order}"]
        osc = PTROscillator(func, samplerate, frequency)
        # pulsar = Pulsar(samplerate, frequency)

        wav = numpy.empty(samplerate)
        # freq = numpy.geomspace(1000, 10000, len(wav))

        for i in range(len(wav)):
            # if pulsar.process() > 0.5:
            #     osc.setPhase(pulsar.phi, pulsar.T)
            # osc.setFrequency(freq[i])
            wav[i] = osc.process()

        snd_dir = Path("snd")
        if not snd_dir.exists():
            snd_dir.mkdir(parents=True)

        soundfile.write(
            str(snd_dir / f"PTRStep{order:02d}.wav"),
            wav,
            samplerate,
            subtype="FLOAT",
        )

render()
