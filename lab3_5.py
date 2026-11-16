# 3.5 Численное интегрирование вариант 12
import json
import math
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
INPUT_FILE = BASE_DIR / "test3_5.json"


def function(x, function_name):
    if function_name == "x_over_x_cubed_plus_8":
        return x / (x**3 + 8)
    raise ValueError("Неизвестная функция")


def antiderivative(x):
    # первообразная для x/(x^3+8)
    if x <= -2:
        raise ValueError("Логарифм от неположительного числа")
    return (-1/6) * math.log(abs(x + 2)) + (1/12) * math.log(x**2 - 2*x + 4) \
           + (math.sqrt(3) / 2) * math.atan((x - 1) / math.sqrt(3))


def exact_integral(a, b):
    # интеграл по формуле Ньютона-Лейбница
    return antiderivative(b) - antiderivative(a)


def integrate(function_name, a, b, h, method):
    n = round((b - a) / h)
    if n <= 0 or abs(a + n * h - b) > 1e-9:
        raise ValueError("Длина отрезка должна делиться на шаг без остатка")
    if method == "simpson" and n % 2 != 0:
        raise ValueError("Для метода Симпсона требуется чётное число разбиений")
    f = lambda x: function(x, function_name)
    if method == "left":
        return h * sum(f(a + i * h) for i in range(n))
    if method == "right":
        return h * sum(f(a + i * h) for i in range(1, n + 1))
    if method == "middle":
        return h * sum(f(a + (i + 0.5) * h) for i in range(n))
    if method == "trapezoid":
        return h * ((f(a) + f(b)) / 2 + sum(f(a + i * h) for i in range(1, n)))
    if method == "simpson":
        total = f(a) + f(b)
        for i in range(1, n):
            total += (4 if i % 2 else 2) * f(a + i * h)
        return h * total / 3
    raise ValueError("Неизвестный метод")


def runge_romberg(F_h1, F_h2, h1, h2, p):
    return F_h2 + (F_h2 - F_h1) / ((h1 / h2)**p - 1)


def main():
    with INPUT_FILE.open(encoding="utf-8") as file:
        config = json.load(file)
    function_name = config["function"]
    a = float(config["a"])
    b = float(config["b"])
    h1 = float(config["h1"])
    h2 = float(config["h2"])
    precision = int(config.get("output_precision", 10))
    if a >= b or not (h1 > h2 > 0):
        raise ValueError("Нужны a<b и h1>h2>0")
    if function_name == "x_over_x_cubed_plus_8" and a <= -2 <= b:
        raise ValueError("Функция имеет разрыв при x=-2")
    
    I_exact = exact_integral(a, b)
    print(f"Интеграл x/(x^3 + 8) на [{a:g}; {b:g}]")
    print(f"Точное значение: {I_exact:.{precision}g}\n")
    print(f"Шаги: h1={h1:g}, h2={h2:g}")
    
    methods = [
        ("left", "Левые прямоугольники", 1),
        ("right", "Правые прямоугольники", 1),
        ("middle", "Средние прямоугольники", 2),
        ("trapezoid", "Трапеции", 2),
        ("simpson", "Симпсон", 4),
    ]
    
    results = []

    for method, name, p in methods:
        F_h1 = integrate(function_name, a, b, h1, method)
        F_h2 = integrate(function_name, a, b, h2, method)
        refined = runge_romberg(F_h1, F_h2, h1, h2, p)
        abs_err = abs(refined - I_exact)
        results.append((name, refined, abs_err))
        
        print(f"\n{name}, порядок p={p}")
        print(f"  I(h1) = {F_h1:.{precision}g}")
        print(f"  I(h2) = {F_h2:.{precision}g}")
        print(f"  Рунге–Ромберг = {refined:.{precision}g}")
        print(f"  Абсолютная погрешность = {abs_err:.{precision}g}")

    best_method, best_val, best_err = min(results, key=lambda x: x[2])
    print(f"\nНаиболее точный метод: {best_method}")
    print(f"  Значение: {best_val:.{precision}g}")
    print(f"  Погрешность: {best_err:.{precision}g}")


if __name__ == "__main__":
    main()