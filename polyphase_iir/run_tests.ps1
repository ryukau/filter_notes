$configs = @(
  @{ order = 4; phase = 2; fc = 0.125 }
  @{ order = 8; phase = 4; fc = 0.03125 }
  @{ order = 12; phase = 4; fc = 0.02 }
  @{ order = 16; phase = 8; fc = 0.015625 }
  @{ order = 16; phase = 32; fc = 0.0078125 }
)

$compilers = @(
  "cl",
  "clang++"
)

foreach ($cc in $compilers) {
  foreach ($cfg in $configs) {
    $id = "$($cfg.order)_$($cfg.phase)"

    python run_tests.py `
      --ref "data/reference_$id.json" `
      --results "data/results_$(($id))_$cc.json" `
      --compiler $cc `
      --samples 48000 `
      --benchmark `
      --design "butterworth,$($cfg.order),$($cfg.fc),$($cfg.phase)"
  }
}
