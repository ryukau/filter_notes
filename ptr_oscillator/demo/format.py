import re
from pathlib import Path

def to_python_code(order, equations):
    indent = " " * 4
    code = f"def PTR{waveform}{order}(phi, T, h):\n"
    code += indent + "n = phi / T\n"
    for idx, eqn in enumerate(equations):
        if idx == 0:
            code += indent + f"if n >= {order - 1:.1f}:\n"
        else:
            code += indent + f"if {idx - 1:.1f} <= n and n < {idx:.1f}:\n"
        eqn = eqn.replace("^", "**")
        code += indent * 2 + f"return {eqn}\n"
    code += indent + "return 0  # Just in case.\n"
    return code

def to_cpp_code(order, equations):
    indent = " " * 2
    typename = "float"
    code = f"{typename} ptr{waveform}{order}({typename} phi, {typename} T, {typename} h) {{\n"
    code += indent + f"{typename} n = phi / T;\n"
    for idx, eqn in enumerate(equations):
        if idx == 0:
            code += indent + f"if (n >= {typename}({order - 1})) "
        else:
            code += indent + f"if (n < {typename}({idx})) "

        eqn = re.sub(
            r"(\w+)\^(\d+)",
            lambda x: (f"{x.group(1)}*" * int(x.group(2)))[:-1],
            eqn,
        )
        eqn = re.sub(
            r"(\d+)",
            lambda x: f"{typename}({x.group(1)})",
            eqn,
        )

        code += f"return {eqn};\n"
    code += indent + "return 0.0; // Just in case.\n}\n"
    return code

filename = Path("maxima_equations/saw")

with open(filename) as eq_file:
    raw = eq_file.read()

raw = raw.replace("\n", "").replace("[", "\n").replace("],", "").replace("]]", "")
indent = "  "
waveform = filename.stem.capitalize()

order = 0
for line in raw.splitlines():
    if re.match(r"^\s*$", line):
        continue
    code = to_python_code(order, line.split(","))
    print(code)
    order += 1
