$compilers = "clang++", "g++", "cl"
$configs = "4_2", "8_4", "12_4", "16_8", "16_32"

foreach ($cc in $compilers) {
  foreach ($cfg in $configs) {
    python plot_results.py `
      --ref "data/reference_$cfg.json" `
      --results "data/results_${cfg}_$cc.json" `
      --save "_${cfg}_$cc"
  }
}
