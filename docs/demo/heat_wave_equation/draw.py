import sys
import math
import random

import OpenGL.GL as gl
from PyQt5.QtCore import (QSize, QTimer, Qt, pyqtSlot, pyqtSignal)
from PyQt5.QtGui import (QPainter, QPainterPath, QBrush, QColor, QPen,
                         QSurfaceFormat, QScreen, QMouseEvent)
from PyQt5.QtWidgets import (QApplication, QVBoxLayout, QOpenGLWidget,
                             QPushButton, QWidget, QHBoxLayout, QDoubleSpinBox,
                             QLabel)

from wave import Wave1D, HeatWave1D, Wave1DExplicit


class Window(QWidget):
    def __init__(self):
        super(Window, self).__init__()

        self.mainLayout = QVBoxLayout()

        self.createGlWidgetLabels()
        self.createGlWidgets()

        self.spinboxC = self.createSlider("c", 4, 0, 100, 0.1, 1)
        self.spinboxDx = self.createSlider("dx", 0.1, 0.001, 1, 0.001, 3)
        self.spinboxDt = self.createSlider("dt", 1 / 60, 0.001, 1, 0.001, 3)
        self.spinboxAlpha = self.createSlider("alpha", 1, 0, 1, 0.001, 3)
        self.spinboxAttenuation = self.createSlider("attenuation", 1, 0.5, 1,
                                                    0.00001, 5)

        self.buttonReset = QPushButton("Reset")
        self.buttonReset.clicked.connect(self.glWidgetImplicit.reset)
        self.mainLayout.addWidget(self.buttonReset)

        self.setLayout(self.mainLayout)

        self.setWindowTitle("1D Wave")

        self.animate()

    def animate(self):
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.glWidgetImplicit.update)
        self.timer.start(1000 / QApplication.primaryScreen().refreshRate())

    def createGlWidgetLabels(self):
        layoutGlLabel = QHBoxLayout()

        label = QLabel("1D Heat Wave")
        label.setAlignment(Qt.AlignCenter)
        layoutGlLabel.addWidget(label)

        self.mainLayout.addLayout(layoutGlLabel)

    def createGlWidgets(self):
        self.layoutGlWidget = QHBoxLayout()
        self.glWidgetImplicit = GLWidget()
        self.layoutGlWidget.addWidget(self.glWidgetImplicit)
        self.mainLayout.addLayout(self.layoutGlWidget)

    def createSlider(self, label, value, minimum, maximum, step, decimals):
        layout = QHBoxLayout()

        label = QLabel(label)
        label.setAlignment(Qt.AlignRight)
        layout.addWidget(label)

        spinbox = QDoubleSpinBox()
        spinbox.setValue(value)
        spinbox.setRange(minimum, maximum)
        spinbox.setSingleStep(step)
        spinbox.setDecimals(decimals)
        layout.addWidget(spinbox)

        self.mainLayout.addLayout(layout)
        spinbox.valueChanged.connect(self.setParameters)

        return spinbox

    def setParameters(self):
        self.setWave1dParameters(self.glWidgetImplicit)

    def setWave1dParameters(self, widget):
        widget.wave1d.setParameters(
            self.spinboxC.value(),
            self.spinboxDx.value(),
            self.spinboxDt.value(),
            self.spinboxAlpha.value(),
            self.spinboxAttenuation.value(),
        )


class GLWidget(QOpenGLWidget):
    mousePressed = pyqtSignal(QMouseEvent)
    mouseMoved = pyqtSignal(QMouseEvent)
    mouseReleased = pyqtSignal(QMouseEvent)

    def __init__(self, parent=None):
        super(GLWidget, self).__init__(parent)
        self.wave1d = HeatWave1D(512, 4, 0.1, 1 / 60, 1, 1)
        self.mousedown = False

    def minimumSizeHint(self):
        return QSize(50, 50)

    def sizeHint(self):
        return QSize(400, 400)

    def move(self):
        self.wave1d.step()

    def paintGL(self):
        self.move()

        size = self.size()
        center = math.floor(size.height() / 2)

        qp = QPainter(self)
        qp.fillRect(0, 0, size.width(), size.height(), Qt.white)

        color = QColor(0x44, 0x88, 0xff)
        qp.setPen(color)
        qp.setBrush(color)

        wav = self.wave1d.value()
        path = QPainterPath()
        path.moveTo(0, center)
        denom = len(wav) - 1
        for i in range(0, len(wav)):
            path.lineTo(
                size.width() * i / denom,
                center + wav[i],
            )
        path.lineTo(size.width() + 1, size.height())
        path.lineTo(0, size.height())
        qp.drawPath(path)

    @pyqtSlot(QMouseEvent)
    def mousePressAction(self, event):
        self.mousedown = True
        self.wave1d.pick(
            event.x() / self.size().width(),
            event.y() - self.size().height() / 2,
        )

    def mousePressEvent(self, event):
        self.mousePressAction(event)
        self.mousePressed.emit(event)

    @pyqtSlot(QMouseEvent)
    def mouseMoveAction(self, event):
        if not self.mousedown:
            return
        self.wave1d.pick(
            event.x() / self.size().width(),
            event.y() - self.size().height() / 2,
        )

    def mouseMoveEvent(self, event):
        self.mouseMoveAction(event)
        self.mouseMoved.emit(event)

    @pyqtSlot(QMouseEvent)
    def mouseReleaseAction(self, event):
        self.mousedown = False
        self.wave1d.pick(
            event.x() / self.size().width(),
            0,
        )

    def mouseReleaseEvent(self, event):
        self.mouseReleaseAction(event)
        self.mouseReleased.emit(event)

    def reset(self):
        self.wave1d.reset()


def antialias(px):
    fmt = QSurfaceFormat()
    fmt.setSamples(px)
    QSurfaceFormat.setDefaultFormat(fmt)


if __name__ == '__main__':
    app = QApplication(sys.argv)

    antialias(10)

    window = Window()
    window.show()
    sys.exit(app.exec_())
