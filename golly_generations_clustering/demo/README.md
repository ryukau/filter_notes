# Clustering of Short Glitch Sound
These scripts make clusters of short (200ms) glitch sounds by using MFCC as feature vectors. Then draw some fancy plots.

# Dependency
- [matplotlib](https://matplotlib.org/)
- [SciPy](https://www.scipy.org/)
- [scikit-learn](https://scikit-learn.org/stable/)
- [PySoundFile](https://pysoundfile.readthedocs.io/en/0.9.0/)
- [python_speech_features](https://python-speech-features.readthedocs.io/en/latest/)

# Usage
Start clustering.

```bash
$ ./start.sh /path/to/data_dir
```

Plot errors for elbow method.

```bash
$ python3 plot_elbow_method.py
```

The repo contains small `test` data sets. `generations` data set can be downloaded from the link below.

- [Golly generations data set](https://drive.google.com/file/d/1wbyGz6bbGULksH3vOchL7za1bPzNSNHs/view?usp=sharing)
