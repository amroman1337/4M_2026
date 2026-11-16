# 3.3 Метод наименьших квадратов вариант 12
import json
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
INPUT_FILE = BASE_DIR / "test3_3.json"


def gauss_solve(matrix, right):
    a = [row[:] for row in matrix]
    b = right[:]
    n = len(b)
    for k in range(n):
        pivot = max(range(k, n), key=lambda i: abs(a[i][k]))
        a[k], a[pivot] = a[pivot], a[k]
        b[k], b[pivot] = b[pivot], b[k]
        if abs(a[k][k]) < 1e-13:
            raise ValueError("Система вырождена, проверьте число и расположение узлов")
        for i in range(k + 1, n):
            multiplier = a[i][k] / a[k][k]
            for j in range(k, n):
                a[i][j] -= multiplier * a[k][j]
            b[i] -= multiplier * b[k]
    solution = [0.0] * n
    for i in range(n - 1, -1, -1):
        solution[i] = (b[i] - sum(a[i][j] * solution[j] for j in range(i + 1, n))) / a[i][i]
    return solution


def fit_polynomial(x, y, degree):    # система нормальных уравнений
    # sum x^(i+j) * a_j = sum y*x^i
    k = degree + 1
    a = [[sum(t ** (i + j) for t in x) for j in range(k)] for i in range(k)]
    b = [sum(value * t**i for t, value in zip(x, y)) for i in range(k)]
    return gauss_solve(a, b)


def polynomial_value(coefficients, x):
    result = 0.0
    for coefficient in reversed(coefficients):
        result = result * x + coefficient
    return result


def polynomial_text(coefficients, precision):
    pieces = []
    for power, value in enumerate(coefficients):
        term = f"{abs(value):.{precision}g}" + ("" if power == 0 else ("*x" if power == 1 else f"*x^{power}"))
        pieces.append(("-" if value < 0 else "+", term))
    return " ".join(
        ("-" if sign == "-" else "") + term if i == 0 else f"{sign} {term}"
        for i, (sign, term) in enumerate(pieces)
    )


def make_plot(x, y, results, filename):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("График пропущен: pip install matplotlib")
        return
    grid = [min(x) + (max(x) - min(x)) * i / 599 for i in range(600)]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.scatter(x, y, color="black", label="Табличные данные")
    for degree, coef in results.items():
        ax.plot(grid, [polynomial_value(coef, t) for t in grid], label=f"Степень {degree}",)
    ax.grid(True)
    ax.legend()
    ax.set(xlabel="x", ylabel="y", title="3.3 Метод наименьших квадратов")
    fig.tight_layout()
    fig.savefig(BASE_DIR / filename, dpi=160)
    plt.close(fig)
    print(f"График создан: {filename}")


def main():
    with INPUT_FILE.open(encoding="utf-8") as file:
        config = json.load(file)
    x = [float(v) for v in config["x"]]
    y = [float(v) for v in config["y"]]
    degrees = [int(v) for v in config.get("degrees", [1, 2, 3])]
    precision = int(config.get("output_precision", 10))
    if len(x) != len(y) or not x or any(d < 0 or d >= len(x) for d in degrees):
        raise ValueError("Проверьте размеры данных и степени полиномов")
    print("Табличные данные:")
    for i, (xi, yi) in enumerate(zip(x, y)):
        print(f"  i={i}, x={xi:g}, y={yi:g}")
    results = {}
    for degree in degrees:
        coefficients = fit_polynomial(x, y, degree)
        results[degree] = coefficients
        deviations = [yi - polynomial_value(coefficients, xi) for xi, yi in zip(x, y)]
        sse = sum(v * v for v in deviations)
        print(f"\nСтепень {degree}")
        print("P(x) =", polynomial_text(coefficients, precision))
        print(f"Сумма квадратов отклонений = {sse:.{precision}g}")
        print(f"Средняя квадратичная ошибка MSE = {sse/len(x):.{precision}g}")
    make_plot(x, y, results, config["plot_file"])


if __name__ == "__main__":
    main()
