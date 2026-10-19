# -*- coding: cp1251 -*-

import sys
import os
import math


def print_matrix(matrix, name="Матрица"):
    print(f"\n{name}:")
    for row in matrix:
        formatted = []
        for val in row:
            if abs(val - round(val)) < 1e-6:
                formatted.append(f"{int(round(val)):>10}")
            else:
                formatted.append(f"{val:>10.6f}")
        print("  [" + " ".join(formatted) + " ]")


def print_vector(vector, name="Вектор"):
    print(f"\n{name}:")
    formatted = []
    for val in vector:
        if abs(val - round(val)) < 1e-6:
            formatted.append(f"{int(round(val)):>10}")
        else:
            formatted.append(f"{val:>10.6f}")
    print("  [" + " ".join(formatted) + " ]")


def read_matrix(filename):
    if not os.path.exists(filename):
        print(f"err: файл {filename} не найден")
        sys.exit(1)

    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    if not lines:
        print("err: файл пустой")
        sys.exit(1)

    try:
        n = int(lines[0])
        if len(lines) < 1 + n:
            print("err: недостаточно данных")
            sys.exit(1)

        A = []
        for i in range(1, n + 1):
            row = list(map(float, lines[i].split()))
            if len(row) != n:
                print(f"err: в строке {i} неверное количество элементов")
                sys.exit(1)
            A.append(row)

        return n, A
    except ValueError:
        print("err: неверный формат данных")
        sys.exit(1)


def max_off_diagonal(n, A):  # максимальный по модулю внедиагональный элемент для критерия остановки
    m = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            if abs(A[i][j]) > m:
                m = abs(A[i][j])
    return m


def jacobi_rotation(n, A, eps, max_iter=1000):
    A = [row[:] for row in A]
    V = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

    iterations = 0
    for _ in range(max_iter):
        max_val = 0.0
        p, q = 0, 1
        for i in range(n):
            for j in range(i + 1, n):
                if abs(A[i][j]) > max_val:
                    max_val = abs(A[i][j])
                    p, q = i, j

        if max_val < eps:
            break

        if abs(A[p][q]) < 1e-15:
            break

        if abs(A[p][p] - A[q][q]) < 1e-15:
            fi = math.pi / 4
        else:
            fi = 0.5 * math.atan(2.0 * A[p][q] / (A[p][p] - A[q][q]))

        c = math.cos(fi)
        s = math.sin(fi)

        A_old = [row[:] for row in A]

        for k in range(n):
            if k == p or k == q:
                continue
            A[k][p] = c * A_old[k][p] + s * A_old[k][q]
            A[k][q] = -s * A_old[k][p] + c * A_old[k][q]
            A[p][k] = A[k][p]
            A[q][k] = A[k][q]

        A[p][p] = c * c * A_old[p][p] + 2.0 * s * c * A_old[p][q] + s * s * A_old[q][q]
        A[q][q] = s * s * A_old[p][p] - 2.0 * s * c * A_old[p][q] + c * c * A_old[q][q]
        A[p][q] = 0.0
        A[q][p] = 0.0

        for k in range(n):
            V_kp = V[k][p]
            V_kq = V[k][q]
            V[k][p] = c * V_kp + s * V_kq
            V[k][q] = -s * V_kp + c * V_kq

        iterations += 1

    eigenvalues = [A[i][i] for i in range(n)]
    return eigenvalues, V, A, iterations


def normalize_signs(n, V):
    for j in range(n):
        for i in range(n):
            if abs(V[i][j]) > 1e-9:
                if V[i][j] < 0:
                    for k in range(n):
                        V[k][j] = -V[k][j]
                break
    return V


def sort_eigenpairs_desc(eigenvalues, V):
    n = len(eigenvalues)
    idx = sorted(range(n), key=lambda i: eigenvalues[i], reverse=True)
    eig_sorted = [eigenvalues[i] for i in idx]
    V_sorted = [[V[i][j] for j in idx] for i in range(n)]
    return eig_sorted, V_sorted


def multiply_matrices(A, B):
    n = len(A)
    m = len(B[0])
    k = len(B)
    C = [[0.0] * m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            for r in range(k):
                C[i][j] += A[i][r] * B[r][j]
    return C


def main():
    eps = 0.001
    filename = "test4.txt"
    n, A = read_matrix(filename)

    print("Задача 1.4: Метод вращений Якоби")
    print(f"Размерность: {n}")
    print_matrix(A, "Исходная матрица A")
    print(f"\n Заданная точность eps = {eps}")

    eigenvalues, V, A_diag, iterations = jacobi_rotation(n, A, eps)
    V = normalize_signs(n, V)
    eigenvalues, V = sort_eigenpairs_desc(eigenvalues, V)

    print(f"\n Собственные значения")
    print(f"  Количество итераций: {iterations}")
    for i in range(n):
        if abs(eigenvalues[i] - round(eigenvalues[i])) < 1e-6:
            print(f"  lambda_{i + 1} = {int(round(eigenvalues[i]))}")
        else:
            print(f"  lambda_{i + 1} = {eigenvalues[i]:.6f}")

    print(f"\n Собственные векторы")
    for j in range(n):
        col = [V[i][j] for i in range(n)]
        print_vector(col, f"h_{j + 1}")
    print_matrix(V, "Матрица собственных векторов (столбцы): ")

    print(f"\n Проверка A * V = V * Lambda")
    Lambda = [[0.0] * n for _ in range(n)]
    for i in range(n):
        Lambda[i][i] = eigenvalues[i]

    AV = multiply_matrices(A, V)
    VLambda = multiply_matrices(V, Lambda)
    print_matrix(AV, "A * V")
    print_matrix(VLambda, "V * Lambda")


if __name__ == "__main__":
    main()