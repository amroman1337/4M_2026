# 3.1  вариант 12
import json
import math
from fractions import Fraction
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
INPUT_FILE = BASE_DIR / "test3_1.json"


def function(x):
    return math.sin(x) + x


def divided_differences(x, y):
    # таблица разделённых разностей
    n = len(x)
    table = [[None] * n for _ in range(n)]
    for i in range(n):
        table[i][0] = y[i]
    for order in range(1, n):
        for i in range(n - order):
            table[i][order] = (table[i + 1][order - 1] - table[i][order - 1]) / (x[i + order] - x[i])
    return table


def lagrange_value(x, y, point):
    result = 0.0
    for i in range(len(x)):
        part = y[i]
        for j in range(len(x)):
            if i != j:
                part *= (point - x[j]) / (x[i] - x[j])
        result += part
    return result


def newton_value(x, coefficients, point):
    result = coefficients[-1]
    for i in range(len(coefficients) - 2, -1, -1):
        result = coefficients[i] + (point - x[i]) * result
    return result


def fraction_label(value):
    fraction = Fraction(*value)
    if fraction == 0:
        return "0"
    top = "" if fraction.numerator == 1 else str(fraction.numerator)
    bottom = "" if fraction.denominator == 1 else f"/{fraction.denominator}"
    return f"{top}π{bottom}"


def polynomial_text(terms):
    pieces = []
    for value, factors in terms:
        sign = "-" if value < 0 else "+"
        text = f"{abs(value):.10g}" + "".join(factors)
        pieces.append((sign, text))
    return " ".join(
        ("-" if sign == "-" else "") + text if i == 0 else f"{sign} {text}"
        for i, (sign, text) in enumerate(pieces)
    )


def print_differences(table, precision):
    print("Таблица разделённых разностей:")
    for i, row in enumerate(table):
        print(f"  i={i}: " + "  ".join("—" if v is None else f"{v:.{precision}g}" for v in row))


def make_plot(x, y, x_star, coefficients, name):
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("График пропущен: pip install matplotlib")
        return
    grid = [x[0] + (x[-1] - x[0]) * i / 599 for i in range(600)]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(grid, [function(t) for t in grid], label="sin(x) + x")
    ax.plot(grid, [lagrange_value(x, y, t) for t in grid], "--", label="Лагранж")
    ax.plot(grid, [newton_value(x, coefficients, t) for t in grid], ":", label="Ньютон")
    ax.scatter(x, y, color="black", label="Узлы")
    ax.scatter([x_star], [lagrange_value(x, y, x_star)], color="red", label="x*")
    ax.set(xlabel="x", ylabel="y", title="задание 3.1")
    ax.grid(True)
    ax.legend()
    fig.tight_layout()
    fig.savefig(BASE_DIR / name, dpi=160)
    plt.close(fig)
    print(f"График: {name}")


def main():
    with INPUT_FILE.open(encoding="utf-8") as file:
        config = json.load(file)
    if config["function"] != "sin_plus_x":
        raise ValueError("Поддерживается только функция sin_plus_x")
    x_star = float(config["x_star"])
    precision = int(config.get("output_precision", 10))
    print(f"f(x) = sin(x) + x; x* = {x_star:g}; точное f(x*) = {function(x_star):.{precision}g}")
    for case in config["cases"]:
        fractions = [Fraction(*p) for p in case["nodes_pi_fractions"]]
        x = [float(frac) * math.pi for frac in fractions]
        if len(x) < 2 or any(x[i] >= x[i + 1] for i in range(len(x) - 1)):
            raise ValueError("Нужно не менее двух возрастающих различных узлов")
        y = [function(t) for t in x]
        table = divided_differences(x, y)
        coeff = table[0]
        labels = [fraction_label(p) for p in case["nodes_pi_fractions"]]
        print(f"\n случай {case['name']}")
        for i, (label, xi, yi) in enumerate(zip(labels, x, y)):
            print(f"  x{i} = {label} = {xi:.{precision}g}; y{i} = {yi:.{precision}g}")
        lag_terms = []
        for i in range(len(x)):
            den = math.prod(x[i] - x[j] for j in range(len(x)) if j != i)
            lag_terms.append((y[i] / den, [f"(x - {labels[j]})" for j in range(len(x)) if j != i]))
        print("L(x) =", polynomial_text(lag_terms))
        print_differences(table, precision)
        new_terms = [(coeff[i], [f"(x - {labels[j]})" for j in range(i)]) for i in range(len(x))]
        print("N(x) =", polynomial_text(new_terms))
        lv = lagrange_value(x, y, x_star)
        nv = newton_value(x, coeff, x_star)
        print(f"L(x*) = {lv:.{precision}g}")
        print(f"N(x*) = {nv:.{precision}g}")
        print(f"|f(x*) - L(x*)| = {abs(function(x_star)-lv):.{precision}g}")
        print(f"|L(x*) - N(x*)| = {abs(lv-nv):.{precision}g}")
        make_plot(x, y, x_star, coeff, case["plot_file"])


if __name__ == "__main__":
    main()
