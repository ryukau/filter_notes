import numpy
import soundfile
from pathlib import Path

# yapf: disable
def PTRSaw0(phi, T, h):
    n = phi / T
    if n >= -1.0:
        return 2*T*n-1
    return 0  # Just in case.

def PTRSaw1(phi, T, h):
    n = phi / T
    if n >= 0.0:
        return 2*T*n-T-1
    if 0.0 <= n and n < 1.0:
        return -2*h*n+2*T*n+2*h-T-1
    return 0  # Just in case.

def PTRSaw2(phi, T, h):
    n = phi / T
    if n >= 1.0:
        return 2*T*n-2*T-1
    if 0.0 <= n and n < 1.0:
        return -h*n**2+2*T*n+2*h-2*T-1
    if 1.0 <= n and n < 2.0:
        return h*n**2-4*h*n+2*T*n+4*h-2*T-1
    return 0  # Just in case.

def PTRSaw3(phi, T, h):
    n = phi / T
    if n >= 2.0:
        return 2*T*n-3*T-1
    if 0.0 <= n and n < 1.0:
        return -(h*n**3)/3+2*T*n+2*h-3*T-1
    if 1.0 <= n and n < 2.0:
        return (2*h*n**3)/3-3*h*n**2+3*h*n+2*T*n+h-3*T-1
    if 2.0 <= n and n < 3.0:
        return -(h*n**3)/3+3*h*n**2-9*h*n+2*T*n+9*h-3*T-1
    return 0  # Just in case.

def PTRSaw4(phi, T, h):
    n = phi / T
    if n >= 3.0:
        return 2*T*n-4*T-1
    if 0.0 <= n and n < 1.0:
        return -(h*n**4)/12+2*T*n+2*h-4*T-1
    if 1.0 <= n and n < 2.0:
        return (h*n**4)/4-(4*h*n**3)/3+2*h*n**2-(4*h*n)/3+2*T*n+(7*h)/3-4*T-1
    if 2.0 <= n and n < 3.0:
        return -(h*n**4)/4+(8*h*n**3)/3-10*h*n**2+(44*h*n)/3+2*T*n-(17*h)/3-4*T-1
    if 3.0 <= n and n < 4.0:
        return (h*n**4)/12-(4*h*n**3)/3+8*h*n**2-(64*h*n)/3+2*T*n+(64*h)/3-4*T-1
    return 0  # Just in case.

def PTRSaw5(phi, T, h):
    n = phi / T
    if n >= 4.0:
        return 2*T*n-5*T-1
    if 0.0 <= n and n < 1.0:
        return -(h*n**5)/60+2*T*n+2*h-5*T-1
    if 1.0 <= n and n < 2.0:
        return (h*n**5)/15-(5*h*n**4)/12+(5*h*n**3)/6-(5*h*n**2)/6+(5*h*n)/12+2*T*n+(23*h)/12-5*T-1
    if 2.0 <= n and n < 3.0:
        return -(h*n**5)/10+(5*h*n**4)/4-(35*h*n**3)/6+(25*h*n**2)/2-(155*h*n)/12+2*T*n+(29*h)/4-5*T-1
    if 3.0 <= n and n < 4.0:
        return (h*n**5)/15-(5*h*n**4)/4+(55*h*n**3)/6-(65*h*n**2)/2+(655*h*n)/12+2*T*n-(133*h)/4-5*T-1
    if 4.0 <= n and n < 5.0:
        return -(h*n**5)/60+(5*h*n**4)/12-(25*h*n**3)/6+(125*h*n**2)/6-(625*h*n)/12+2*T*n+(625*h)/12-5*T-1
    return 0  # Just in case.

def PTRSaw6(phi, T, h):
    n = phi / T
    if n >= 5.0:
        return 2*T*n-6*T-1
    if 0.0 <= n and n < 1.0:
        return -(h*n**6)/360+2*T*n+2*h-6*T-1
    if 1.0 <= n and n < 2.0:
        return (h*n**6)/72-(h*n**5)/10+(h*n**4)/4-(h*n**3)/3+(h*n**2)/4-(h*n)/10+2*T*n+(121*h)/60-6*T-1
    if 2.0 <= n and n < 3.0:
        return -(h*n**6)/36+(2*h*n**5)/5-(9*h*n**4)/4+(19*h*n**3)/3-(39*h*n**2)/4+(79*h*n)/10+2*T*n-(13*h)/20-6*T-1
    if 3.0 <= n and n < 4.0:
        return (h*n**6)/36-(3*h*n**5)/5+(21*h*n**4)/4-(71*h*n**3)/3+(231*h*n**2)/4-(731*h*n)/10+2*T*n+(797*h)/20-6*T-1
    if 4.0 <= n and n < 5.0:
        return -(h*n**6)/72+(2*h*n**5)/5-(19*h*n**4)/4+(89*h*n**3)/3-(409*h*n**2)/4+(1829*h*n)/10+2*T*n-(7849*h)/60-6*T-1
    if 5.0 <= n and n < 6.0:
        return (h*n**6)/360-(h*n**5)/10+(3*h*n**4)/2-12*h*n**3+54*h*n**2-(648*h*n)/5+2*T*n+(648*h)/5-6*T-1
    return 0  # Just in case.

def PTRSaw7(phi, T, h):
    n = phi / T
    if n >= 6.0:
        return 2*T*n-7*T-1
    if 0.0 <= n and n < 1.0:
        return -(h*n**7)/2520+2*T*n+2*h-7*T-1
    if 1.0 <= n and n < 2.0:
        return (h*n**7)/420-(7*h*n**6)/360+(7*h*n**5)/120-(7*h*n**4)/72+(7*h*n**3)/72-(7*h*n**2)/120+(7*h*n)/360+2*T*n+(719*h)/360-7*T-1
    if 2.0 <= n and n < 3.0:
        return -(h*n**7)/168+(7*h*n**6)/72-(77*h*n**5)/120+(161*h*n**4)/72-(329*h*n**3)/72+(133*h*n**2)/24-(1337*h*n)/360+2*T*n+(1103*h)/360-7*T-1
    if 3.0 <= n and n < 4.0:
        return (h*n**7)/126-(7*h*n**6)/36+(119*h*n**5)/60-(98*h*n**4)/9+(1253*h*n**3)/36-(196*h*n**2)/3+(12089*h*n)/180+2*T*n-(1229*h)/45-7*T-1
    if 4.0 <= n and n < 5.0:
        return -(h*n**7)/168+(7*h*n**6)/36-(161*h*n**5)/60+(182*h*n**4)/9-(3227*h*n**3)/36+(700*h*n**2)/3-(59591*h*n)/180+2*T*n+(9011*h)/45-7*T-1
    if 5.0 <= n and n < 6.0:
        return (h*n**7)/420-(7*h*n**6)/72+(203*h*n**5)/120-(1169*h*n**4)/72+(6671*h*n**3)/72-(7525*h*n**2)/24+(208943*h*n)/360+2*T*n-(162287*h)/360-7*T-1
    if 6.0 <= n and n < 7.0:
        return -(h*n**7)/2520+(7*h*n**6)/360-(49*h*n**5)/120+(343*h*n**4)/72-(2401*h*n**3)/72+(16807*h*n**2)/120-(117649*h*n)/360+2*T*n+(117649*h)/360-7*T-1
    return 0  # Just in case.

def PTRSaw8(phi, T, h):
    n = phi / T
    if n >= 7.0:
        return 2*T*n-8*T-1
    if 0.0 <= n and n < 1.0:
        return -(h*n**8)/20160+2*T*n+2*h-8*T-1
    if 1.0 <= n and n < 2.0:
        return (h*n**8)/2880-(h*n**7)/315+(h*n**6)/90-(h*n**5)/45+(h*n**4)/36-(h*n**3)/45+(h*n**2)/90-(h*n)/315+2*T*n+(5041*h)/2520-8*T-1
    if 2.0 <= n and n < 3.0:
        return -(h*n**8)/960+(2*h*n**7)/105-(13*h*n**6)/90+(3*h*n**5)/5-(55*h*n**4)/36+(37*h*n**3)/15-(223*h*n**2)/90+(149*h*n)/105+2*T*n+(829*h)/504-8*T-1
    if 3.0 <= n and n < 4.0:
        return (h*n**8)/576-(h*n**7)/21+(5*h*n**6)/9-(18*h*n**5)/5+(128*h*n**4)/9-(106*h*n**3)/3+(488*h*n**2)/9-(4954*h*n)/105+2*T*n+(6259*h)/315-8*T-1
    if 4.0 <= n and n < 5.0:
        return -(h*n**8)/576+(4*h*n**7)/63-h*n**6+(398*h*n**5)/45-48*h*n**4+(1474*h*n**3)/9-344*h*n**2+(128498*h*n)/315+2*T*n-(7269*h)/35-8*T-1
    if 5.0 <= n and n < 6.0:
        return (h*n**8)/960-(h*n**7)/21+(17*h*n**6)/18-(53*h*n**5)/5+(2647*h*n**4)/36-(967*h*n**3)/3+(15683*h*n**2)/18-(139459*h*n)/105+2*T*n+(2211007*h)/2520-8*T-1
    if 6.0 <= n and n < 7.0:
        return -(h*n**8)/2880+(2*h*n**7)/105-(41*h*n**6)/90+(31*h*n**5)/5-(1889*h*n**4)/36+(4237*h*n**3)/15-(84881*h*n**2)/90+(187133*h*n)/105+2*T*n-(3667649*h)/2520-8*T-1
    if 7.0 <= n and n < 8.0:
        return (h*n**8)/20160-(h*n**7)/315+(4*h*n**6)/45-(64*h*n**5)/45+(128*h*n**4)/9-(4096*h*n**3)/45+(16384*h*n**2)/45-(262144*h*n)/315+2*T*n+(262144*h)/315-8*T-1
    return 0  # Just in case.

def PTRSaw9(phi, T, h):
    n = phi / T
    if n >= 8.0:
        return 2*T*n-9*T-1
    if 0.0 <= n and n < 1.0:
        return -(h*n**9)/181440+2*T*n+2*h-9*T-1
    if 1.0 <= n and n < 2.0:
        return (h*n**9)/22680-(h*n**8)/2240+(h*n**7)/560-(h*n**6)/240+(h*n**5)/160-(h*n**4)/160+(h*n**3)/240-(h*n**2)/560+(h*n)/2240+2*T*n+(40319*h)/20160-9*T-1
    if 2.0 <= n and n < 3.0:
        return -(h*n**9)/6480+(h*n**8)/320-(3*h*n**7)/112+(31*h*n**6)/240-(63*h*n**5)/160+(127*h*n**4)/160-(17*h*n**3)/16+(73*h*n**2)/80-(1023*h*n)/2240+2*T*n+(42367*h)/20160-9*T-1
    if 3.0 <= n and n < 4.0:
        return (h*n**9)/3240-(3*h*n**8)/320+(69*h*n**7)/560-(221*h*n**6)/240+(693*h*n**5)/160-(2141*h*n**4)/160+(2183*h*n**3)/80-(2843*h*n**2)/80+(60213*h*n)/2240+2*T*n-(141341*h)/20160-9*T-1
    if 4.0 <= n and n < 5.0:
        return -(h*n**9)/2592+(h*n**8)/64-(31*h*n**7)/112+(45*h*n**6)/16-(2891*h*n**5)/160+(2439*h*n**4)/32-(10159*h*n**3)/48+(5985*h*n**2)/16-(857291*h*n)/2240+2*T*n+(78415*h)/448-9*T-1
    if 5.0 <= n and n < 6.0:
        return (h*n**9)/3240-(h*n**8)/64+(39*h*n**7)/112-(215*h*n**6)/48+(5859*h*n**5)/160-(6311*h*n**4)/32+(11197*h*n**3)/16-(25265*h*n**2)/16+(4611459*h*n)/2240+2*T*n-(4763015*h)/4032-9*T-1
    if 6.0 <= n and n < 7.0:
        return -(h*n**9)/6480+(3*h*n**8)/320-(141*h*n**7)/560+(941*h*n**6)/240-(6237*h*n**5)/160+(41021*h*n**4)/160-(89167*h*n**3)/80+(246923*h*n**2)/80-(11064957*h*n)/2240+2*T*n+(70243421*h)/20160-9*T-1
    if 7.0 <= n and n < 8.0:
        return (h*n**9)/22680-(h*n**8)/320+(11*h*n**7)/112-(431*h*n**6)/240+(3367*h*n**5)/160-(26207*h*n**4)/160+(40619*h*n**3)/48-(223673*h*n**2)/80+(11994247*h*n)/2240+2*T*n-(91171007*h)/20160-9*T-1
    if 8.0 <= n and n < 9.0:
        return -(h*n**9)/181440+(h*n**8)/2240-(9*h*n**7)/560+(27*h*n**6)/80-(729*h*n**5)/160+(6561*h*n**4)/160-(19683*h*n**3)/80+(531441*h*n**2)/560-(4782969*h*n)/2240+2*T*n+(4782969*h)/2240-9*T-1
    return 0  # Just in case.

def PTRSaw10(phi, T, h):
    n = phi / T
    if n >= 9.0:
        return 2*T*n-10*T-1
    if 0.0 <= n and n < 1.0:
        return -(h*n**10)/1814400+2*T*n+2*h-10*T-1
    if 1.0 <= n and n < 2.0:
        return (h*n**10)/201600-(h*n**9)/18144+(h*n**8)/4032-(h*n**7)/1512+(h*n**6)/864-(h*n**5)/720+(h*n**4)/864-(h*n**3)/1512+(h*n**2)/4032-(h*n)/18144+2*T*n+(362881*h)/181440-10*T-1
    if 2.0 <= n and n < 3.0:
        return -(h*n**10)/50400+(h*n**9)/2268-(17*h*n**8)/4032+(5*h*n**7)/216-(71*h*n**6)/864+(143*h*n**5)/720-(287*h*n**4)/864+(575*h*n**3)/1512-(1151*h*n**2)/4032+(329*h*n)/2592+2*T*n+(358273*h)/181440-10*T-1
    if 3.0 <= n and n < 4.0:
        return (h*n**10)/21600-(h*n**9)/648+(13*h*n**8)/576-(289*h*n**7)/1512+(901*h*n**6)/864-(2773*h*n**5)/720+(8461*h*n**4)/864-(3667*h*n**3)/216+(11083*h*n**2)/576-(233893*h*n)/18144+2*T*n+(1066861*h)/181440-10*T-1
    if 4.0 <= n and n < 5.0:
        return -(h*n**10)/14400+(h*n**9)/324-(35*h*n**8)/576+(1055*h*n**7)/1512-(4475*h*n**6)/864+(18731*h*n**5)/720-(77555*h*n**4)/864+(45485*h*n**3)/216-(185525*h*n**2)/576+(5271131*h*n)/18144+2*T*n-(4190647*h)/36288-10*T-1
    if 5.0 <= n and n < 6.0:
        return (h*n**10)/14400-(5*h*n**9)/1296+(55*h*n**8)/576-(2095*h*n**7)/1512+(11275*h*n**6)/864-(60019*h*n**5)/720+(316195*h*n**4)/864-(235765*h*n**3)/216+(1220725*h*n**2)/576-(43947619*h*n)/18144+2*T*n+(45028103*h)/36288-10*T-1
    if 6.0 <= n and n < 7.0:
        return -(h*n**10)/21600+(h*n**9)/324-(53*h*n**8)/576+(2441*h*n**7)/1512-(15941*h*n**6)/864+(103277*h*n**5)/720-(663581*h*n**4)/864+(604043*h*n**3)/216-(3818123*h*n**2)/576+(167683997*h*n)/18144+2*T*n-(1044649181*h)/181440-10*T-1
    if 7.0 <= n and n < 8.0:
        return (h*n**10)/50400-(h*n**9)/648+(31*h*n**8)/576-(1675*h*n**7)/1512+(12871*h*n**6)/864-(98407*h*n**5)/720+(748207*h*n**4)/864-(807745*h*n**3)/216+(6064393*h*n**2)/576-(316559287*h*n)/18144+2*T*n+(2345053807*h)/181440-10*T-1
    if 8.0 <= n and n < 9.0:
        return -(h*n**10)/201600+(h*n**9)/2268-(71*h*n**8)/4032+(629*h*n**7)/1512-(5561*h*n**6)/864+(49049*h*n**5)/720-(431441*h*n**4)/864+(3782969*h*n**3)/1512-(33046721*h*n**2)/4032+(287420489*h*n)/18144+2*T*n-(2486784401*h)/181440-10*T-1
    if 9.0 <= n and n < 10.0:
        return (h*n**10)/1814400-(h*n**9)/18144+(5*h*n**8)/2016-(25*h*n**7)/378+(125*h*n**6)/108-(125*h*n**5)/9+(3125*h*n**4)/27-(125000*h*n**3)/189+(156250*h*n**2)/63-(3125000*h*n)/567+2*T*n+(3125000*h)/567-10*T-1
    return 0  # Just in case.
# yapf: enable

class PTROscillator:
    def __init__(self, ptr_func, samplerate, frequency):
        self.ptr_func = ptr_func
        self.samplerate = samplerate
        self.setFrequency(frequency)
        self.phi = 0.0
        self.h = 1.0

    def setFrequency(self, frequency):
        self.T = frequency / self.samplerate

    def setPhase(self, phi, T):
        self.phi = phi
        ratio = self.T / T
        self.h = ratio - numpy.floor(ratio)

    def process(self):
        self.phi += self.T
        if self.phi >= 1.0:
            self.h = 1.0
            self.phi -= 1.0
        return self.ptr_func(self.phi, self.T, self.h)

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

def renderPureSaw():
    def render(directory, prefix, phase, T):
        for order in range(11):
            func = globals()[f"PTRSaw{order}"]
            ptr = [func(phi, T, 1) for phi in phase]
            soundfile.write(str(directory / f"{prefix}{order:02d}.wav"), ptr, samplerate)

    samplerate = 44100
    frequency = 1000
    T = frequency / samplerate
    sync_ratio = 4 / 3

    phase = numpy.linspace(0, frequency, samplerate) % 1.0

    snd_dir = Path("snd")
    if not snd_dir.exists():
        snd_dir.mkdir(parents=True)

    render(snd_dir, "PTRSaw", phase, T)

def renderSyncSaw():
    samplerate = 44100
    frequency = 1000
    sync_ratio = 3 / 2

    for order in range(11):
        func = globals()[f"PTRSaw{order}"]
        osc = PTROscillator(func, samplerate, sync_ratio * frequency)
        pulsar = Pulsar(samplerate, frequency)

        wav = numpy.empty(samplerate)
        freq = numpy.geomspace(1000, 10000, len(wav))

        for i in range(len(wav)):
            if pulsar.process() > 0.5:
                osc.setPhase(pulsar.phi, pulsar.T)
            osc.setFrequency(freq[i])
            wav[i] = osc.process()

        snd_dir = Path("snd")
        prefix = f"Sync{sync_ratio}"
        if not snd_dir.exists():
            snd_dir.mkdir(parents=True)

        soundfile.write(
            str(snd_dir / f"PTRSawSync{order:02d}.wav"),
            wav,
            samplerate,
            subtype="FLOAT",
        )

renderPureSaw()
renderSyncSaw()
