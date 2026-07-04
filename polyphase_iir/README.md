To run tests:

```ps1
run_tests.ps1 # or `run_tests.sh`
plot_results.ps1
```

## `run_tests.py`
Run tests on C++ implementations. `python run_tests.py --help` to list the options.

Invocation examples:

```ps1
python run_tests.py --compiler cl --benchmark --design "butterworth,8,0.3125,4"
python run_tests.py --compiler clang++ --benchmark --design "butterworth,2,0.125,4"
python run_tests.py --compiler g++ --benchmark --design "elliptic,8,0.1,60,0.2,4"
```

## `plot_results.py`
Plot the results of `run_tests.py`. `python plot_results.py --help` to list the options.

## Other Files
To design a polyphase IIR, `design.py` and `signal_mp.py` are required.
