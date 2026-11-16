# 3.4 Численное дифференцирование вариант 12
import json
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
INPUT_FILE = BASE_DIR / "test3_4.json"


def three_point_derivatives(x, y, point):
    if len(x) < 3 or len(x) != len(y) or len(set(x)) != len(x):
        raise ValueError("Нужны минимум три различных узла")
    indexes = sorted(range(len(x)), key=lambda i: (abs(x[i] - point), i))[:3]
    indexes.sort(key=lambda i: x[i])
    xx = [x[i] for i in indexes]
    yy = [y[i] for i in indexes]
    first = 0.0
    second = 0.0
    for i in range(3):
        j, k = [t for t in range(3) if t != i]
        den = (xx[i] - xx[j]) * (xx[i] - xx[k])
        first += yy[i] * (2 * point - xx[j] - xx[k]) / den
        second += 2 * yy[i] / den
    return first, second, xx, yy


def adjacent_derivatives(x, y, point):
    if point not in x:
        return None
    i = x.index(point)
    if i == 0 or i == len(x) - 1:
        return None
    left = (y[i] - y[i - 1]) / (x[i] - x[i - 1])
    right = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
    return left, right


def make_plot(x, y, point, filename):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("График пропущен: pip install matplotlib")
        return
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(x, y, "o-", label="Табличные значения")
    ax.axvline(point, color="red", linestyle="--", label="x*")
    ax.grid(True)
    ax.legend()
    ax.set(xlabel="x", ylabel="y", title="3.4 Численное дифференцирование")
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
    if any(x[i] >= x[i + 1] for i in range(len(x) - 1)):
        raise ValueError("x должен строго возрастать")
    print("Исходная таблица:")
    for i, (xi, yi) in enumerate(zip(x, y)):
        print(f"  i={i}, x={xi:g}, y={yi:g}")
    first, second, xx, yy = three_point_derivatives(x, y, point)
    print(f"\nТочка x* = {point:g}")
    print("Три узла для квадратичного полинома:")
    for xi, yi in zip(xx, yy):
        print(f"  ({xi:g}; {yi:g})")
    print(f"y'(x*) ≈ {first:.{precision}g}")
    print(f"y''(x*) ≈ {second:.{precision}g}")
    slopes = adjacent_derivatives(x, y, point)
    if slopes is not None:
        print(f"Наклон слева ≈ {slopes[0]:.{precision}g}")
        print(f"Наклон справа ≈ {slopes[1]:.{precision}g}")
    make_plot(x, y, point, config["plot_file"])


if __name__ == "__main__":
    main()
