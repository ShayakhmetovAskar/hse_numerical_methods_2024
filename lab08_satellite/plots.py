import pandas as pd
import matplotlib.pyplot as plt

# Чтение данных из CSV файла
data = pd.read_csv("satellite_orbit.csv")

# Построение графика траектории спутника
plt.figure(figsize=(10, 10))
plt.plot(data['PositionX'], data['PositionY'], label='Trajectory')
plt.xlabel('X Position (km)')
plt.ylabel('Y Position (km)')
plt.title('Satellite Orbit Trajectory')
plt.legend()
plt.grid(True)
plt.axis('equal')
plt.show()

# Построение графика изменения высоты
data['Radius'] = (data['PositionX']**2 + data['PositionY']**2 + data['PositionZ']**2)**0.5
plt.figure(figsize=(10, 6))
plt.plot(data['Step'], data['Radius'], label='Altitude')
plt.xlabel('Time Step (s)')
plt.ylabel('Altitude (km)')
plt.title('Satellite Altitude Over Time')
plt.legend()
plt.grid(True)
plt.show()
