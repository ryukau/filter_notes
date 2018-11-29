import sys
import math
import random
import pathlib

from PyQt5.QtCore import (QSize, QPointF, Qt)
from PyQt5.QtGui import (QGuiApplication, QPainter, QPainterPath, QBrush,
                         QColor, QPen, QImage, QFont)

from wave import Wave1D, HeatWave1D, Wave1DExplicit


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
    qp.drawText(4, font.pointSize() + 4, "α={:0.3}".format(wave1d.alpha))

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

    width = 320
    height = 180

    grid = 4
    length = grid * grid
    denom = length - 1
    wave1d = [
        HeatWave1D(
            width,
            4,
            0.1 * width / 512,
            1 / 60,
            index / denom,
            1,
        ) for index in range(length)
    ]

    images = [
        QImage(width, height, QImage.Format_RGB888) for _ in range(length)
    ]
    composed_image = QImage(width * grid, height * grid, QImage.Format_RGB888)

    font = QFont("Dejavu Sans Mono", 12)

    pathlib.Path('img').mkdir(parents=True, exist_ok=True)

    for wave in wave1d:
        wave.pick(0.5, height / 3)
    num_frame = 1800
    for time in range(num_frame):
        if time == num_frame / 6:
            for wave in wave1d:
                wave.pick(0.5, 0)
        for i in range(length):
            wave1d[i].step()
            draw(images[i], wave1d[i], font)
        compose(composed_image, images, grid, width, height)
        composed_image.save("img/out{:04d}.png".format(time))
