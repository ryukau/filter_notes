import sys
import math
import random
import pathlib

from PyQt5.QtCore import (QSize, QPointF, Qt)
from PyQt5.QtGui import (QGuiApplication, QPainter, QPainterPath, QBrush,
                         QColor, QPen, QImage, QFont, QFontMetrics)

from wave import *


def draw(image, wave1d, font):
    size = image.size()
    centerY = math.floor(size.height() / 2)

    qp = QPainter(image)
    qp.setRenderHints(QPainter.Antialiasing)
    qp.fillRect(0, 0, size.width(), size.height(), Qt.white)

    color = QColor(0x44, 0x88, 0xff)
    qp.setPen(color)
    qp.setBrush(color)

    wav = wave1d.value()
    path = QPainterPath()
    path.moveTo(0, centerY)
    denom = len(wav) - 1
    for i in range(0, len(wav)):
        path.lineTo(
            size.width() * i / denom,
            centerY + wav[i],
        )
    path.lineTo(size.width() + 1, size.height())
    path.lineTo(0, size.height())
    qp.drawPath(path)

    if wave1d.pickY != 0:
        pen = QPen(QColor(0xe0, 0x20, 0x10))
        pen.setWidthF(2.0)
        pen.setCapStyle(Qt.RoundCap)
        qp.setPen(pen)
        w = 8
        h = 16
        centerX = size.width() / 2 - 0.5
        depth = centerY + wave1d.pickY
        qp.drawLine(QPointF(centerX, centerY - h), QPointF(centerX, centerY))
        qp.drawLine(
            QPointF(centerX - w / 2, centerY - h / 4),
            QPointF(centerX, centerY),
        )
        qp.drawLine(
            QPointF(centerX + w / 2, centerY - h / 4),
            QPointF(centerX, centerY),
        )

    qp.setPen(Qt.black)
    qp.setFont(font)
    line_space = QFontMetrics(font).lineSpacing()
    qp.drawText(4, 1 * line_space, "  α={:0.3}".format(wave1d.alpha))
    qp.drawText(4, 2 * line_space, "  β={:0.3}".format(wave1d.beta))
    qp.drawText(4, 3 * line_space, "τ_ε={:0.8}".format(wave1d.tau_epsilon))

    qp.setPen(QColor(0x30, 0x30, 0x30))
    qp.setBrush(Qt.transparent)
    qp.drawRect(0, 0, size.width(), size.height())


def compose(target, images, grid, gridWidth, gridHeight):
    size = target.size()
    qp = QPainter(target)
    for x in range(grid):
        for y in range(grid):
            qp.drawImage(
                gridWidth * x,
                gridHeight * y,
                images[x + grid * y],
            )


if __name__ == '__main__':
    app = QGuiApplication(sys.argv)

    grid = 4
    length = grid * grid
    denom = length - 1

    width = int(1280 / grid)
    height = int(720 / grid)

    wave1d = [
        ZenerWave1D(
            width,
            16,
            0.1 * width / 512,
            1 / 60,
            1,
            1,
            1,
            2 * i / denom,
            2 * i / denom,
        ) for i in range(length)
    ]

    pathlib.Path('img').mkdir(parents=True, exist_ok=True)

    images = [
        QImage(width, height, QImage.Format_RGB888) for _ in range(length)
    ]
    composed_image = QImage(width * grid, height * grid, QImage.Format_RGB888)

    font = QFont("Dejavu Sans Mono", 12)

    current_frame = 0
    jMax = 5
    for i in range(3):
        for j in range(jMax):
            for wave in wave1d:
                wave.reset()
                wave.tau_epsilon = 0.1**(i + 1)
                wave.beta = j * 2 / (jMax - 1)
                wave.pick(0.5, height / 6)
            num_frame = 1200
            for time in range(num_frame):
                if time == num_frame / 6:
                    for wave in wave1d:
                        wave.pick(0.5, 0)
                for k in range(length):
                    wave1d[k].step()
                    draw(images[k], wave1d[k], font)
                compose(composed_image, images, grid, width, height)
                composed_image.save("img/out{:08d}.png".format(current_frame))
                sys.stdout.write(
                    "\rScene: {:2}.{:2}, Frame: {:4} / {:4}".format(
                        i, j, time, num_frame))
                current_frame += 1
    print("\n")
