import numpy as np
import matplotlib.pyplot as plt
import skia
import scipy.signal as signal

paint_background = skia.Paint(
    Color=skia.ColorWHITE,
    Style=skia.Paint.kFill_Style,
)
paint_sig = skia.Paint(
    AntiAlias=True,
    Style=skia.Paint.kStroke_Style,
    Color=skia.ColorBLUE,
    StrokeWidth=2,
)
paint_discarded_sig = skia.Paint(
    AntiAlias=True,
    Style=skia.Paint.kStroke_Style,
    Color=skia.ColorGRAY,
    StrokeWidth=1,
)
paint_split = skia.Paint(
    AntiAlias=True,
    PathEffect=skia.DashPathEffect.Make([5.0, 5.0, 5.0, 5.0], 0.0),
    Style=skia.Paint.kStroke_Style,
    Color=skia.Color4f(0x88888888),
    StrokeWidth=2,
)

paint_text = skia.Paint(Style=skia.Paint.kFill_Style, Color=skia.ColorBLACK)
style = skia.FontStyle(
    skia.FontStyle.kBold_Weight,
    skia.FontStyle.kCondensed_Width,
    skia.FontStyle.kUpright_Slant,
)
font = skia.Font(skia.Typeface("Dejavu Sans", style), 16)

def signalToPath(sig, x, y, width, gain):
    p = skia.Path()
    values = -gain * sig + y
    xscale = width / len(values)
    p.moveTo(x, values[0])
    for idx, value in enumerate(values):
        p.lineTo(x + xscale * idx, value)
    return p

def plotOverlapAddSplit(blockSize, nBlock, length, source, fir, blocks):
    width = 520
    height = 600
    surface = skia.Surface(width, height)

    with surface as cv:
        cv.drawRect(skia.Rect(0, 0, width, height), paint_background)

        # FIR.
        cv.drawTextBlob(skia.TextBlob("FIR", font), 10, 100, paint_text)
        path = signalToPath(16 * fir, 100, 100, 100, 32)
        cv.drawPath(path, paint_sig)

        # Source.
        cv.drawTextBlob(skia.TextBlob("Source", font), 10, 200, paint_text)
        path = signalToPath(source, 100, 200, 300, 32)
        cv.drawPath(path, paint_sig)

        # Blocks.
        for idx in range(nBlock):
            left = (idx + 1) * 100
            yCenter = (idx + 3) * 100
            cv.drawTextBlob(skia.TextBlob(f"Block {idx}", font), 10, yCenter, paint_text)

            path = signalToPath(blocks[idx], left, yCenter, 200, 32)
            cv.drawPath(path, paint_sig)

        # Vertical split.
        for idx in range(nBlock + 1):
            left = (idx + 1) * 100
            path = skia.Path()
            path.moveTo(left, 0)
            path.lineTo(left, height)
            cv.drawPath(path, paint_split)

    surface.makeImageSnapshot().save("img/overlap_add_split.png", skia.kPNG)

def plotOverlapAddMerge(blockSize, nBlock, length, merged, fir, blocks):
    width = 620
    height = 600
    surface = skia.Surface(width, height)

    with surface as cv:
        cv.drawRect(skia.Rect(0, 0, width, height), paint_background)

        # FIR.
        cv.drawTextBlob(skia.TextBlob("FIR", font), 10, 100, paint_text)
        path = signalToPath(16 * fir, 150, 100, 100, 32)
        cv.drawPath(path, paint_sig)

        # Signal blocks.
        for idx in range(nBlock):
            left = (idx + 1.5) * 100
            yCenter = (idx + 2) * 100
            cv.drawTextBlob(skia.TextBlob(f"FIR * Block {idx}", font), 10, yCenter,
                            paint_text)
            path = signalToPath(blocks[idx], left, yCenter, 200, 32)
            cv.drawPath(path, paint_sig)

        # Operators.
        for idx in range(nBlock - 1):
            cv.drawTextBlob(skia.TextBlob("+", font), 10, (idx + 2.5) * 100, paint_text)
        cv.drawTextBlob(skia.TextBlob("=", font), 10, (nBlock + 1.5) * 100, paint_text)

        # Vertical split.
        for idx in range(nBlock + 1):
            left = (idx + 1.5) * 100
            path = skia.Path()
            path.moveTo(left, 0)
            path.lineTo(left, height)
            cv.drawPath(path, paint_split)

        # Merged.
        cv.drawTextBlob(skia.TextBlob("Output", font), 10, 500, paint_text)
        path = signalToPath(merged, 150, 500, 400, 32)
        cv.drawPath(path, paint_sig)

    surface.makeImageSnapshot().save("img/overlap_add_merge.png", skia.kPNG)

def plotOverlapAdd():
    blockSize = 1024
    nBlock = 3
    length = nBlock * blockSize
    source = signal.sawtooth(np.linspace(0, 2 * np.pi * 10, length))
    fir = signal.firwin(blockSize, 960, fs=48000)

    spc_fir = np.fft.rfft(fir, 2 * blockSize)

    source_blocks = []
    filtered_blocks = []
    merged = np.zeros((nBlock + 1) * blockSize)
    for i in range(nBlock):
        splitted = np.hstack([
            source[i * blockSize:(i + 1) * blockSize],
            np.zeros(blockSize),
        ])
        source_blocks.append(splitted)

        filtered = np.fft.irfft(np.fft.rfft(splitted) * spc_fir)
        filtered_blocks.append(filtered)

        merged[i * blockSize:(i + 2) * blockSize] += filtered

    plotOverlapAddSplit(blockSize, nBlock, length, source, fir, source_blocks)
    plotOverlapAddMerge(blockSize, nBlock, length, merged, fir, filtered_blocks)

def plotOverlapSaveSplit(blockSize, nBlock, length, source, fir, blocks):
    width = 520
    height = 600
    surface = skia.Surface(width, height)

    with surface as cv:
        cv.drawRect(skia.Rect(0, 0, width, height), paint_background)

        # FIR.
        cv.drawTextBlob(skia.TextBlob("FIR", font), 10, 100, paint_text)
        path = signalToPath(16 * fir, 100, 100, 100, 32)
        cv.drawPath(path, paint_sig)

        # Source.
        cv.drawTextBlob(skia.TextBlob("Source", font), 10, 200, paint_text)
        # path = signalToPath(source[:blockSize], 100, 200, 100, 32)
        # cv.drawPath(path, paint_discarded_sig)
        path = signalToPath(source[blockSize:-blockSize], 200, 200, 200, 32)
        cv.drawPath(path, paint_sig)
        # path = signalToPath(source[-blockSize:], 400, 200, 100, 32)
        # cv.drawPath(path, paint_discarded_sig)

        # Blocks.
        for idx in range(nBlock):
            left = (idx + 1) * 100
            yCenter = (idx + 3) * 100
            cv.drawTextBlob(skia.TextBlob(f"Block {idx}", font), 10, yCenter, paint_text)
            path = signalToPath(blocks[idx], left, yCenter, 200, 32)
            cv.drawPath(path, paint_sig)

        # Vertical split.
        for idx in range(nBlock):
            left = (idx + 2) * 100
            path = skia.Path()
            path.moveTo(left, 0)
            path.lineTo(left, height)
            cv.drawPath(path, paint_split)

    surface.makeImageSnapshot().save("img/overlap_save_split.png", skia.kPNG)

def plotOverlapSaveMerge(blockSize, nBlock, length, merged, fir, blocks):
    width = 620
    height = 600
    surface = skia.Surface(width, height)

    with surface as cv:
        cv.drawRect(skia.Rect(0, 0, width, height), paint_background)

        # FIR.
        cv.drawTextBlob(skia.TextBlob("FIR", font), 10, 100, paint_text)
        path = signalToPath(16 * fir, 150, 100, 100, 32)
        cv.drawPath(path, paint_sig)

        # Signal blocks.
        for idx in range(nBlock):
            left = (idx + 1.5) * 100
            yCenter = (idx + 2) * 100
            cv.drawTextBlob(skia.TextBlob(f"FIR * Block {idx}", font), 10, yCenter,
                            paint_text)

            # Discarded part.
            path = signalToPath(blocks[idx][:blockSize], left, yCenter, 100, 32)
            cv.drawPath(path, paint_discarded_sig)

            # Valid part.
            path = signalToPath(blocks[idx][blockSize:], left + 100, yCenter, 100, 32)
            cv.drawPath(path, paint_sig)

        # Operators.
        for idx in range(nBlock - 1):
            cv.drawTextBlob(skia.TextBlob("+", font), 10, (idx + 2.5) * 100, paint_text)
        cv.drawTextBlob(skia.TextBlob("=", font), 10, (nBlock + 1.5) * 100, paint_text)

        # Vertical split.
        for idx in range(nBlock):
            left = (idx + 2.5) * 100
            path = skia.Path()
            path.moveTo(left, 0)
            path.lineTo(left, height)
            cv.drawPath(path, paint_split)

        # Merged.
        cv.drawTextBlob(skia.TextBlob("Output", font), 10, 500, paint_text)
        path = signalToPath(merged[:blockSize], 150, 500, 100, 32)
        cv.drawPath(path, paint_discarded_sig)
        path = signalToPath(merged[blockSize:], 250, 500, 300, 32)
        cv.drawPath(path, paint_sig)

    surface.makeImageSnapshot().save("img/overlap_save_merge.png", skia.kPNG)

def plotOverlapSave():
    blockSize = 1024
    nBlock = 3
    length = nBlock * blockSize
    source = np.hstack([
        np.zeros(blockSize),
        signal.sawtooth(np.linspace(0, 2 * np.pi * 7, length - blockSize)),
        np.zeros(blockSize),
    ])
    fir = signal.firwin(blockSize, 960, fs=48000)

    spc_fir = np.fft.rfft(fir, 2 * blockSize)

    source_blocks = []
    filtered_blocks = []
    merged = np.zeros((nBlock + 1) * blockSize)
    for i in range(nBlock):
        splitted = source[i * blockSize:(i + 2) * blockSize]
        source_blocks.append(splitted)

        filtered = np.fft.irfft(np.fft.rfft(splitted) * spc_fir)
        filtered_blocks.append(filtered)

        merged[(i + 1) * blockSize:(i + 2) * blockSize] += filtered[blockSize:]

    plotOverlapSaveSplit(blockSize, nBlock, length, source, fir, source_blocks)
    plotOverlapSaveMerge(blockSize, nBlock, length, merged, fir, filtered_blocks)

if __name__ == "__main__":
    plotOverlapAdd()
    plotOverlapSave()
