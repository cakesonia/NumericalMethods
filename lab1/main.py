import numpy.linalg as la
import numpy as np
import math
from scipy import reshape


class Matrix:
    def __init__(self):
        k = 19
        p = 21
        s = 0.02 * k
        b = 0.02 * p
        self.n = 4
        self.A = [[8.3, 2.62 + s, 4.1, 1.9],
                  [3.92, 8.45, 7.78 - s, 2.46],
                  [3.77, 7.21 + s, 8.04, 2.28],
                  [2.21, 3.65 - s, 1.69, 6.69]]
        self.B = [-10.65 + b, 12.21, 15.45 - b, -8.35]
        self.V = range(self.n * self.n)
        self.V = reshape(self.V, (self.n, self.n))
        self.P = range(self.n * self.n)
        self.P = reshape(self.P, (self.n, self.n))
        self.C = range(self.n)
        self.C = reshape(self.C, (self.n, 1))
        self.inx = 0

    def straight_move(self):
        for i in range(0, self.n - 1):
            self.inx = i
        for i, j in range(0, self.n - 1):
            self.V[i][j] = self.A[i][j]
            self.P[i] = self.B[i]
        for k in range(0, self.n - 1):
            max = math.fabs(self.V[k][k])
            w = k
            for i in range(k, self.n - 1):
                if max < math.fabs(self.V[k][i]):
                    max = math.fabs(self.V[k][i])
                    w = i
            z = self.inx[k]
            self.inx[k] = self.inx[w]
            self.inx[w] = z
            for d in range(0, self.n - 1):
                if d < k:
                    value = self.C[d][k]
                    self.C[d][k] = self.C[d][w]
                    self.C[d][w] = value
                else:
                    value = self.V[d][k]
                    self.V[d][k] = self.V[d][w]
                    self.V[d][w] = value
            Y = range(self.n)
            Y = reshape(Y, (self.n, 1))
            Y[k] = self.P[k] / self.V[k][k]
            for i in range(k, self.n - 1):
                self.P[i] = self.P[i] - (self.V[i][k] * Y[k])
                for j in range(k, self.n - 1):
                    self.C[k][j] = self.V[k][j] / self.V[k][k]
                    self.V[j][j] = self.V[i][j] - (self.V[i][k] * self.C[k][j])
            return la.solve(self.A, self.B)

    def reverse_move(self, Y):
        X = range(self.n)
        X = reshape(X, (self.n, 1))
        X[self.n - 1] = float(Y[self.n - 1])
        for i in range(self.n - 2, 0, -1):
            sum = 0
            for j in range(self.n - 1, i - 1, -1):
                sum += self.C[i][j] * X[j]
            X[i] = float((Y[i] - sum) / self.C[i][i])
        return X

    def solve():
        A = [[8.3, 3, 4.1, 1.9],
             [3.92, 8.45, 7.4, 2.46],
             [3.77, 7.59, 8.04, 2.28],
             [2.21, 3.27, 1.69, 6.69]]
        B = [-10.23, 12.21, 11.25, -8.35]
        X = la.solve(A, B)
        print(X)
        C = np.dot(A, X)
        if np.array_equal(C, B):
            pass
        return C

    def ordering(self, X):
        for i in range(0, self.n - 1):
            if self.inx != i:
                z = self.inx[i]
                value = X[i]
                X[i] = X[z]
                X[z] = value
                self.inx[i] = self.inx[z]
                self.inx[z] = z
            return X


if __name__ == "__main__":
    print(Matrix.solve())
