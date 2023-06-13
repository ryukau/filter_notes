`resonance.py` searches max resonance of `FeedbackEMALowpass`. It takes long time to finish, and outputs `table_linear.json` or `table_log.json`. Following is meaning of values in `*.json`:

- `normalizedFreq`: Cutoff frequency in \[rad/2π\].
- `resonance`: Maximum resonance value of `FeedbackEMALowpass` at `normalizedFreq` on same index.

`plot.py` shows cutoff-resonance plot from `table.json`. Rename `table_linear.json` or `table_log.json`.
