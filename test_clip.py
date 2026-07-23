import sys
from PySide6.QtWidgets import QApplication, QLabel, QVBoxLayout, QWidget, QPushButton
from PySide6.QtGui import QPixmap, Qt

app = QApplication(sys.argv)
w = QWidget()
l = QVBoxLayout(w)

lbl = QLabel()
px = QPixmap("/Users/jzuluaga/dev/MontuPython/montu_gui/assets/montu-python-logo-complete.png").scaledToWidth(180, Qt.TransformationMode.SmoothTransformation)
lbl.setPixmap(px)
lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
lbl.setMinimumSize(px.size())
l.addWidget(lbl)

for i in range(15):
    l.addWidget(QPushButton(f"Button {i}"))

w.show()
print("Min height of window:", w.minimumHeight())
sys.exit(0)
