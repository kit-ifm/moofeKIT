import numpy as np
import matplotlib as mlp
import matplotlib.pyplot as plt
import tikzplotlib

# matlab: x = [1:1:10]
x = np.arange(0,11,1,float)
y = np.ones_like(x)
#x = np.arange(0, 10, 1)  # Start, Stop, Step
#y = np.sin(x)
# Einen x-y-Plot erstellen:
#plt.plot([1,2,3,4], [1,4,9,16], 'bo-')
plt.plot(x, y, 'b-')

# Label für die x-Achse vergeben:
plt.xlabel('time')

# Label für die y-Achse vergeben:
plt.ylabel('energy')

# Achsen-Bereiche manuell festlegen
# Syntax: plt.axis([xmin, xmax, ymin, ymax])
#plt.axis([0, 5, 0, 20])

# Ein gepunktetes Diagramm-Gitter einblenden:
#plt.grid(True)

plt.savefig("energy.png",dpi=300)
plt.savefig("energy.eps",dpi=300)
tikzplotlib.save("mytikz.tex")

# Diagramm anzeigen:
plt.show()

