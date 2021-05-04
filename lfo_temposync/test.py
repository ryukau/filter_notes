import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import soundfile
import gc

def calcDistance(diff):
    if diff > 0.5:  # ラップアラウンドあり。target は後ろにいる。
        return diff - 1
    if diff < -0.5:  # ラップアラウンドあり。target は前にいる。
        return 1 + diff
    return diff  # ラップアラウンドなし。

def generateSin(fs, freq, duration, nFrame):
    phase = np.linspace(0, 2 * np.pi * freq * duration, nFrame)
    return np.sin(phase)

def naiveSync(beatsElapsed, syncInterval):
    phi = beatsElapsed / syncInterval
    return phi - np.floor(phi)

def plotSpectrogram(sig, samplerate, name):
    freq, time, spectre = signal.spectrogram(
        sig,
        samplerate,
        # ("tukey", 0.33),
        "blackmanharris",
        nperseg=2048,
        noverlap=2048 - 2048 // 4,
    )

    absed = np.abs(spectre)
    spectre = 20 * np.log10(absed / np.max(absed))
    spectre = np.where(spectre < -200, -200, spectre)

    plt.figure(figsize=(8, 4), tight_layout=True)
    plt.pcolormesh(time, freq, spectre, cmap="magma", shading="gouraud")
    plt.yscale("log")
    plt.ylim((20, samplerate / 2))
    plt.ylabel("Frequency [Hz]")
    plt.xlabel("Time [s]")
    plt.savefig("img/" + name + ".png")
    plt.close("all")

    gc.collect()  # Maybe out of memory on 32bit CPython.

def testSyncIntervalChange():
    fs = 48000
    duration = 2
    nFrame = int(fs * duration)
    nFrame += nFrame % 2
    time = np.linspace(0, duration, nFrame)

    tempo = np.full(nFrame, 120)
    sync = np.hstack([
        np.full(nFrame // 2, 6 / 5),
        np.full(nFrame // 2, 2),
    ])
    beats = (tempo / (fs * 60)).cumsum()

    plt.figure(figsize=(6, 3), tight_layout=True)
    plt.plot(time, sync, color="black")
    plt.title("Sync Interval")
    plt.xlabel("Time [s]")
    plt.ylabel("Sync Interval [beat]")
    plt.ylim((0.95, 2.05))
    plt.grid()
    plt.savefig("img/SyncChangeSyncIntervalStep.svg")
    plt.clf()

    lfo = np.empty(nFrame)
    for idx in range(len(sync)):
        lfo[idx] = naiveSync(beats[idx], sync[idx])

    plt.figure(figsize=(6, 3), tight_layout=True)
    plt.plot(time, lfo, color="black")
    plt.title("LFO Phase")
    plt.xlabel("Time [s]")
    plt.ylabel("Normalized Phase")
    plt.grid()
    plt.savefig("img/SyncChangeLfoPhaseStep.svg")
    plt.clf()

    sm = 0.5 * np.sin(2 * np.pi * lfo) + 0.5

    plt.figure(figsize=(6, 3), tight_layout=True)
    plt.plot(time, sm, color="black")
    plt.title("Sin LFO")
    plt.xlabel("Time [s]")
    plt.ylabel("Amplitude")
    plt.grid()
    plt.savefig("img/SyncChangeAmModulatorStep.svg")
    plt.clf()

    sc = generateSin(fs, 1000, duration, nFrame)

    sig = sm * sc

    plotSpectrogram(sig, fs, "NaiveSyncChange")

    soundfile.write("snd/NaiveSyncChange.wav", sig, fs, subtype="FLOAT")

def testTempoChange():
    fs = 48000
    duration = 2
    nFrame = int(fs * duration)
    nFrame += nFrame % 2
    time = np.linspace(0, duration, nFrame)

    tempo = np.hstack([
        np.full(nFrame // 2, 40),
        np.full(nFrame // 2, 120),
    ])
    sync = np.full(nFrame, 1)
    beats = (tempo / (fs * 60)).cumsum()

    plt.figure(figsize=(6, 3), tight_layout=True)
    plt.plot(time, sync, color="black")
    plt.title("Sync Interval")
    plt.xlabel("Time [s]")
    plt.ylabel("Sync Interval [beat]")
    plt.ylim((0.95, 2.05))
    plt.grid()
    plt.savefig("img/TempoChangeSyncIntervalStep.svg")
    plt.clf()

    lfo = np.empty(nFrame)
    for idx in range(len(sync)):
        lfo[idx] = naiveSync(beats[idx], sync[idx])

    plt.figure(figsize=(6, 3), tight_layout=True)
    plt.plot(time, lfo, color="black")
    plt.title("LFO Phase")
    plt.xlabel("Time [s]")
    plt.ylabel("Normalized Phase")
    plt.grid()
    plt.savefig("img/TempoChangeLfoPhaseStep.svg")
    plt.clf()

    sm = 0.5 * np.sin(2 * np.pi * lfo) + 0.5

    plt.figure(figsize=(6, 3), tight_layout=True)
    plt.plot(time, sm, color="black")
    plt.title("Sin LFO")
    plt.xlabel("Time [s]")
    plt.ylabel("Amplitude")
    plt.grid()
    plt.savefig("img/TempoChangeAmModulatorStep.svg")
    plt.clf()

    sc = generateSin(fs, 1000, duration, nFrame)

    sig = sm * sc

    plotSpectrogram(sig, fs, "NaiveTempoChange")

    soundfile.write("snd/NaiveTempoChange.wav", sig, fs, subtype="FLOAT")

# testSyncIntervalChange()
testTempoChange()
