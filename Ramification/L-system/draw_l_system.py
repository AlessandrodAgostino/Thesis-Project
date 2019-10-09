from l_system import L_System
from branch import Branch

ls = L_System()
ls.start()
ls.multiple_iterations(8)
ls.draw(c = 'r', circle = 'smart')
ls.savefig('RamificationSmartCircles.png', dpi=1000)
