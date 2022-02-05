import pandas as pd
import requests
import io
import pyBounds

# import data: p_t, q_t
url = "https://raw.githubusercontent.com/JMSLab/PyBounds/870170743f155b3a9da680ec64e0d5079e266601/PyBounds/roberts_schlenker_2013.csv?token=GHSAT0AAAAAABQYAQ7N7L2PIVTKPDHHG7WWYQEWC6A"
download = requests.get(url).content

data = pd.read_csv(io.StringIO(download.decode('utf-8')))
df = pd.DataFrame(data)
p_t = df.iloc[:, 0]
q_t = df.iloc[:, 1]

if __name__ == '__main__':
    bounds_data = pyBounds.Bounds(p_t, q_t)
    print(bounds_data.delta_epsilon(5))