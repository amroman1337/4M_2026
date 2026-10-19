# -*- coding: utf_8 -*-

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


def identity_matrix(n):
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]


def copy_matrix(M):
    return [row[:] for row in M]


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


def vector_norm(v):
    return math.sqrt(sum(x * x for x in v))


def householder_matrix(n, v):
    v_norm_sq = sum(x * x for x in v)
    if v_norm_sq < 1e-30:
        return identity_matrix(n)
    H = identity_matrix(n)
    for i in range(n):
        for j in range(n):
            H[i][j] -= 2.0 * v[i] * v[j] / v_norm_sq
    return H


def qr_decomposition(A):
    n = len(A)
    R = copy_matrix(A)
    Q = identity_matrix(n)

    for k in range(n - 1):
        b = [R[i][k] for i in range(k, n)]
        norm_b = vector_norm(b)

        if norm_b < 1e-15:
            continue

        sign_b1 = 1.0 if b[0] >= 0 else -1.0
        v_full = [0.0] * n
        v_full[k] = b[0] + sign_b1 * norm_b
        for i in range(k + 1, n):
            v_full[i] = b[i - k]

        v = [v_full[i] for i in range(k, n)]

        H_k = householder_matrix(n - k, v)
        H = identity_matrix(n)
        for i in range(k, n):
            for j in range(k, n):
                H[i][j] = H_k[i - k][j - k]

        R = multiply_matrices(H, R)
        Q = multiply_matrices(Q, H)

    return Q, R


def qr_algorithm(n, A, eps=1e-10, max_iter=2000):
    A_k = copy_matrix(A)
    eigenvalues = [0.0] * n

    for iteration in range(max_iter):
        Q, R = qr_decomposition(A_k)
        A_new = multiply_matrices(R, Q)
        diff = max(abs(A_new[i][j] - A_k[i][j]) for i in range(n) for j in range(n))
        A_k = A_new

        if diff < eps:
            break

    for i in range(n):
        eigenvalues[i] = A_k[i][i]

    return eigenvalues, A_k, iteration + 1


def main():
    eps = 0.01
    filename = "test5.txt"
    n, A = read_matrix(filename)

    print("Задача 1.5: QR-разложение и QR-алгоритм")
    print(f"Размерность: {n}")
    print_matrix(A, "Матрица A")
    Q, R = qr_decomposition(A)

    print_matrix(Q, "Матрица Q")
    print_matrix(R, "Матрица R")

    QR = multiply_matrices(Q, R)
    print_matrix(QR, "Результат Q * R:")
    eigenvalues, A_final, iters = qr_algorithm(n, A, eps)

    print(f"\n Собственные значения: {iters} итераций)")
    for i in range(n):
        if abs(eigenvalues[i] - round(eigenvalues[i])) < 1e-6:
            print(f"   lambda_{i + 1} = {int(round(eigenvalues[i]))}")
        else:
            print(f"   lambda_{i + 1} = {eigenvalues[i]:.6f}")

    print_matrix(A_final, "Матрица после QR-итераций")
    return n, A


if __name__ == "__main__":
    n, A = main()
