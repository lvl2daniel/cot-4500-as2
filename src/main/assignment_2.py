import numpy as np
np.set_printoptions(floatmode='fixed')


def neville_interpolation(x_vals, y_vals, target_x):
    n = len(x_vals)
    table = np.zeros((n, n))
    
    for i in range(n):
        table[i][0] = y_vals[i]
    
    for j in range(1, n):
        for i in range(j, n):
            term1 = (target_x - x_vals[i - j]) * table[i][j - 1]
            term2 = (target_x - x_vals[i]) * table[i - 1][j - 1]
            table[i][j] = (term1 - term2) / (x_vals[i] - x_vals[i - j])
    
    return table[n - 1][n - 1]

def newton_divided_differences(x_vals, y_vals):
    n = len(x_vals)
    table = np.zeros((n, n))
    
    for i in range(n):
        table[i][0] = y_vals[i]
    
    for j in range(1, n):
        for i in range(n - j):
            table[i][j] = (table[i + 1][j - 1] - table[i][j - 1]) / (x_vals[i + j] - x_vals[i])
    
    return table

def newton_interpolation(x_vals, divided_diffs, target_x):
    n = len(x_vals)
    result = divided_diffs[0][0]
    t = 1.0
    
    for i in range(n-1):
        t *= (target_x - x_vals[i])
        result += divided_diffs[0][i+1] * t
    
    return result

def hermite_interpolation(x_vals, y_vals, derivatives):
    n = len(x_vals)
    m = 2 * n
    table = np.zeros((m, m))
    
    for i in range(n):
        table[2 * i][0] = x_vals[i]
        table[2 * i + 1][0] = x_vals[i]
        table[2 * i][1] = y_vals[i]
        table[2 * i + 1][1] = y_vals[i]
        table[2 * i + 1][2] = derivatives[i]
    
    for j in range(2, m):
        for i in range(j, m):
            if table[i][j] == 0:
                table[i][j] = (table[i][j - 1] - table[i - 1][j - 1]) / (table[i][0] - table[i - j + 1][0])
    
    return table

def cubic_spline_interpolation(x_vals, y_vals):
    n = len(x_vals)
    h_vals = np.diff(x_vals)
    A = np.zeros((n, n))
    rhs = np.zeros(n)
    
    A[0][0] = A[n - 1][n - 1] = 1
    
    for i in range(1, n - 1):
        A[i][i - 1] = h_vals[i - 1]
        A[i][i] = 2 * (h_vals[i - 1] + h_vals[i])
        A[i][i + 1] = h_vals[i]
        rhs[i] = 3 * ((y_vals[i + 1] - y_vals[i]) / h_vals[i] - (y_vals[i] - y_vals[i - 1]) / h_vals[i - 1])
    
    spline_coeffs = np.linalg.solve(A, rhs)
    
    return A, rhs, spline_coeffs

def main():
    x_vals = [3.6, 3.8, 3.9]
    y_vals = [1.675, 1.436, 1.318]
    target_x = 3.7
    
    nev_result = neville_interpolation(x_vals, y_vals, target_x)
    print("Question 1:\n", nev_result, "\n")
    
    x1_vals = [7.2, 7.4, 7.5, 7.6]
    y1_vals = [23.5492, 25.3913, 26.8224, 27.4589]
    
    newton_table = newton_divided_differences(x1_vals, y1_vals)
    print("Question 2:\n")
    print(f'{newton_table[0][1]}\n{newton_table[0][2]}\n{newton_table[0][3]}\n')
    print("\n")
    
    newton_result = newton_interpolation(x1_vals, newton_table, 7.3)
    print("Question 3:\n",newton_result, "\n")
    
    x3_vals = [3.6, 3.8, 3.9]
    f_vals = [1.675, 1.436, 1.318]
    f_derivatives = [-1.195, -1.188, -1.182]
    
    hermite_result = hermite_interpolation(x3_vals, f_vals, f_derivatives)
    print("Question 4:\n")
    for row in hermite_result:
        print(" ".join(f"{num:.7e}" for num in row))
    print("\n")
    
    x_spline = [2, 5, 8, 10]
    y_spline = [3, 5, 7, 9]
    

    A, b, spline_coeffs = cubic_spline_interpolation(x_spline, y_spline)
    print("Question 5:\n")
    for row in A:
        print(row)
    print("\n", b)
    print("\n", spline_coeffs)
    
if __name__ == "__main__":
    main()
