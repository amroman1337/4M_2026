# 3.2 кубический сплайн вариант 12
import json
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
INPUT_FILE = BASE_DIR / "test3_2.json"


def tridiagonal_solve(lower, diagonal, upper, right):
    n = len(right)
    if n == 0:
        return []
    p = [0.0] * n
    q = [0.0] * n
    p[0] = upper[0] / diagonal[0] if n > 1 else 0.0
    q[0] = right[0] / diagonal[0]
    for i in range(1, n):
        denominator = diagonal[i] - lower[i] * p[i - 1]
        p[i] = upper[i] / denominator if i < n - 1 else 0.0
        q[i] = (right[i] - lower[i] * q[i - 1]) / denominator
    answer = [0.0] * n
    answer[-1] = q[-1]
    for i in range(n - 2, -1, -1):
        answer[i] = q[i] - p[i] * answer[i + 1]
    return answer


def spline_coefficients(x, y):
    n = len(x)
    if n < 2 or len(y) != n or any(x[i] >= x[i + 1] for i in range(n - 1)):
        raise ValueError("массивы x/y не одинаковые, нужны минимум два возрастающих узла")
    h = [x[i + 1] - x[i] for i in range(n - 1)]
    count = n - 2
    lower = [0.0] * count
    diagonal = [0.0] * count
    upper = [0.0] * count
    right = [0.0] * count
    for r in range(count):
        i = r + 1
        lower[r] = h[i - 1] if r > 0 else 0.0
        diagonal[r] = 2 * (h[i - 1] + h[i])
        upper[r] = h[i] if r < count - 1 else 0.0
        right[r] = 3 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1])
    c = [0.0] + tridiagonal_solve(lower, diagonal, upper, right) + [0.0]
    a = []
    b = []
    d = []
    for i in range(n - 1):
        a.append(y[i])
        b.append((y[i + 1] - y[i]) / h[i] - h[i] * (2 * c[i] + c[i + 1]) / 3)
        d.append((c[i + 1] - c[i]) / (3 * h[i]))
    return a, b, c[:-1], d


def spline_value(x, coefficients, point):
    if not x[0] <= point <= x[-1]:
        raise ValueError("Точка должна лежать внутри диапазона узлов")
    a, b, c, d = coefficients
    i = len(x) - 2
    for j in range(len(x) - 1):
        if point <= x[j + 1]:
            i = j
            break
    t = point - x[i]
    return a[i] + t * (b[i] + t * (c[i] + t * d[i]))


def make_plot(x, y, coefficients, point, filename):
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("График пропущен: pip install matplotlib")
        return
    grid = [x[0] + (x[-1] - x[0]) * i / 599 for i in range(600)]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(grid, [spline_value(x, coefficients, t) for t in grid], label="Кубический сплайн")
    ax.scatter(x, y, color="black", label="Исходные узлы")
    ax.scatter([point], [spline_value(x, coefficients, point)], color="red", label="x*")
    ax.grid(True)
    ax.legend()
    ax.set(xlabel="x", ylabel="y", title="3.2. Естественный кубический сплайн")
    fig.tight_layout()
    fig.savefig(BASE_DIR / filename, dpi=160)
    plt.close(fig)
    print(f"График: {filename}")


def main():
    with INPUT_FILE.open(encoding="utf-8") as file:
        config = json.load(file)
    x = [float(v) for v in config["x"]]
    y = [float(v) for v in config["y"]]
    point = float(config["x_star"])
    precision = int(config.get("output_precision", 10))
    coeff = spline_coefficients(x, y)
    a, b, c, d = coeff
    print("Исходные данные:")
    for i, (xi, yi) in enumerate(zip(x, y)):
        print(f"  x{i}={xi:g}, f{i}={yi:.{precision}g}")
    print("\nУсловия: S''(x0) = S''(xn)=0")
    print("Форма: S_i(x)=a_i+b_i*(x-x_(i-1))+c_i*(x-x_(i-1))²+d_i*(x-x_(i-1))³")
    print(f"{'i':>3} {'интервал':>16} {'a_i':>14} {'b_i':>14} {'c_i':>14} {'d_i':>14}")
    for j in range(len(a)):
        print(
            f"{j+1:>3} {f'[{x[j]:g}; {x[j+1]:g}]':>16} "
            f"{a[j]:>14.{precision}g} {b[j]:>14.{precision}g} "
            f"{c[j]:>14.{precision}g} {d[j]:>14.{precision}g}"
        )
    print(f"\nS({point:g}) = {spline_value(x,coeff,point):.{precision}g}")
    max_error = max(abs(spline_value(x, coeff, xi) - yi) for xi, yi in zip(x, y))
    print(f"Проверка прохождения через узлы, max ошибка = {max_error:.3g}")
    make_plot(x, y, coeff, point, config["plot_file"])


if __name__ == "__main__":
    main()
