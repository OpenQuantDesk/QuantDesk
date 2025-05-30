#pragma once

#include "types.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>

namespace math::core {

class Statistics {
public:
    template<typename Iterator>
    static Real mean(Iterator first, Iterator last) {
        if (first == last) return 0.0;
        return std::accumulate(first, last, Real{0}) / std::distance(first, last);
    }
    
    template<typename Iterator>
    static Real variance(Iterator first, Iterator last, Real mean_val) {
        if (first == last) return 0.0;
        Real sum = 0.0;
        Size count = 0;
        for (auto it = first; it != last; ++it) {
            Real diff = *it - mean_val;
            sum += diff * diff;
            ++count;
        }
        return count > 1 ? sum / (count - 1) : 0.0;
    }
    
    template<typename Iterator>
    static Real standardDeviation(Iterator first, Iterator last, Real mean_val) {
        return std::sqrt(variance(first, last, mean_val));
    }
    
    template<typename Iterator>
    static Real percentile(Iterator first, Iterator last, Real p) {
        if (first == last) return 0.0;
        Vector sorted(first, last);
        std::sort(sorted.begin(), sorted.end());
        
        Real index = p * (sorted.size() - 1);
        Size lower = static_cast<Size>(std::floor(index));
        Size upper = static_cast<Size>(std::ceil(index));
        
        if (lower == upper) return sorted[lower];
        
        Real weight = index - lower;
        return sorted[lower] * (1.0 - weight) + sorted[upper] * weight;
    }
    
    template<typename Iterator>
    static Real skewness(Iterator first, Iterator last, Real mean_val, Real std_val) {
        if (first == last || std_val == 0.0) return 0.0;
        Real sum = 0.0;
        Size count = 0;
        for (auto it = first; it != last; ++it) {
            Real normalized = (*it - mean_val) / std_val;
            sum += normalized * normalized * normalized;
            ++count;
        }
        return count > 0 ? sum / count : 0.0;
    }
    
    template<typename Iterator>
    static Real kurtosis(Iterator first, Iterator last, Real mean_val, Real std_val) {
        if (first == last || std_val == 0.0) return 0.0;
        Real sum = 0.0;
        Size count = 0;
        for (auto it = first; it != last; ++it) {
            Real normalized = (*it - mean_val) / std_val;
            Real squared = normalized * normalized;
            sum += squared * squared;
            ++count;
        }
        return count > 0 ? (sum / count) - 3.0 : 0.0;
    }
};

class NumericalMethods {
public:
    template<typename Function>
    static Real bisection(Function f, Real a, Real b, Real tolerance = TOLERANCE, Integer maxIter = MAX_ITERATIONS) {
        Real fa = f(a);
        Real fb = f(b);
        
        if (fa * fb > 0) return std::numeric_limits<Real>::quiet_NaN();
        
        for (Integer iter = 0; iter < maxIter; ++iter) {
            Real c = (a + b) * 0.5;
            Real fc = f(c);
            
            if (std::abs(fc) < tolerance || std::abs(b - a) < tolerance) {
                return c;
            }
            
            if (fa * fc < 0) {
                b = c;
                fb = fc;
            } else {
                a = c;
                fa = fc;
            }
        }
        
        return (a + b) * 0.5;
    }
    
    template<typename Function>
    static Real newtonRaphson(Function f, Function df, Real x0, Real tolerance = TOLERANCE, Integer maxIter = MAX_ITERATIONS) {
        Real x = x0;
        
        for (Integer iter = 0; iter < maxIter; ++iter) {
            Real fx = f(x);
            Real dfx = df(x);
            
            if (std::abs(fx) < tolerance) return x;
            if (std::abs(dfx) < EPSILON) break;
            
            Real newX = x - fx / dfx;
            if (std::abs(newX - x) < tolerance) return newX;
            
            x = newX;
        }
        
        return x;
    }
    
    template<typename Function>
    static Real brent(Function f, Real a, Real b, Real tolerance = TOLERANCE, Integer maxIter = MAX_ITERATIONS) {
        Real fa = f(a);
        Real fb = f(b);
        
        if (fa * fb > 0) return std::numeric_limits<Real>::quiet_NaN();
        
        Real c = a;
        Real fc = fa;
        Real d = b - a;
        Real e = d;
        
        for (Integer iter = 0; iter < maxIter; ++iter) {
            if (std::abs(fc) < std::abs(fb)) {
                a = b; b = c; c = a;
                fa = fb; fb = fc; fc = fa;
            }
            
            Real tol = 2.0 * EPSILON * std::abs(b) + tolerance;
            Real m = 0.5 * (c - b);
            
            if (std::abs(m) <= tol || std::abs(fb) < tolerance) {
                return b;
            }
            
            if (std::abs(e) >= tol && std::abs(fa) > std::abs(fb)) {
                Real s, p, q;
                if (a == c) {
                    s = fb / fa;
                    p = 2.0 * m * s;
                    q = 1.0 - s;
                } else {
                    q = fa / fc;
                    Real r = fb / fc;
                    s = fb / fa;
                    p = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0));
                    q = (q - 1.0) * (r - 1.0) * (s - 1.0);
                }
                
                if (p > 0) q = -q;
                else p = -p;
                
                if (2.0 * p < std::min(3.0 * m * q - std::abs(tol * q), std::abs(e * q))) {
                    e = d;
                    d = p / q;
                } else {
                    d = m;
                    e = d;
                }
            } else {
                d = m;
                e = d;
            }
            
            a = b;
            fa = fb;
            
            if (std::abs(d) > tol) {
                b += d;
            } else {
                b += (m > 0) ? tol : -tol;
            }
            
            fb = f(b);
            
            if (fb * fc > 0) {
                c = a;
                fc = fa;
                d = b - a;
                e = d;
            }
        }
        
        return b;
    }
    
    template<typename Function>
    static Real trapezoidalRule(Function f, Real a, Real b, Integer n = 1000) {
        Real h = (b - a) / n;
        Real sum = 0.5 * (f(a) + f(b));
        
        for (Integer i = 1; i < n; ++i) {
            sum += f(a + i * h);
        }
        
        return sum * h;
    }
    
    template<typename Function>
    static Real simpsonsRule(Function f, Real a, Real b, Integer n = 1000) {
        if (n % 2 == 1) ++n;
        
        Real h = (b - a) / n;
        Real sum = f(a) + f(b);
        
        for (Integer i = 1; i < n; i += 2) {
            sum += 4.0 * f(a + i * h);
        }
        
        for (Integer i = 2; i < n; i += 2) {
            sum += 2.0 * f(a + i * h);
        }
        
        return sum * h / 3.0;
    }
};

class Interpolation {
public:
    static Real linear(Real x, Real x0, Real x1, Real y0, Real y1) {
        if (std::abs(x1 - x0) < EPSILON) return y0;
        return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
    }
    
    static Real bilinear(Real x, Real y, Real x0, Real x1, Real y0, Real y1,
                        Real z00, Real z01, Real z10, Real z11) {
        Real t = (x - x0) / (x1 - x0);
        Real u = (y - y0) / (y1 - y0);
        
        return (1.0 - t) * (1.0 - u) * z00 +
               t * (1.0 - u) * z10 +
               (1.0 - t) * u * z01 +
               t * u * z11;
    }
    
    static Vector cubicSpline(const Vector& x, const Vector& y, const Vector& xi) {
        Size n = x.size();
        if (n < 2 || y.size() != n) return {};
        
        Vector h(n - 1);
        for (Size i = 0; i < n - 1; ++i) {
            h[i] = x[i + 1] - x[i];
        }
        
        Vector alpha(n - 1);
        for (Size i = 1; i < n - 1; ++i) {
            alpha[i] = 3.0 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);
        }
        
        Vector l(n), mu(n), z(n);
        l[0] = 1.0;
        
        for (Size i = 1; i < n - 1; ++i) {
            l[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
        }
        
        l[n - 1] = 1.0;
        z[n - 1] = 0.0;
        
        Vector c(n), b(n - 1), d(n - 1);
        for (Integer j = n - 2; j >= 0; --j) {
            c[j] = z[j] - mu[j] * c[j + 1];
            b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
            d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
        }
        
        Vector result(xi.size());
        for (Size k = 0; k < xi.size(); ++k) {
            Real xval = xi[k];
            Size j = 0;
            
            while (j < n - 1 && x[j + 1] < xval) ++j;
            if (j >= n - 1) j = n - 2;
            
            Real dx = xval - x[j];
            result[k] = y[j] + b[j] * dx + c[j] * dx * dx + d[j] * dx * dx * dx;
        }
        
        return result;
    }
};

class MatrixOperations {
public:
    static Matrix multiply(const Matrix& A, const Matrix& B) {
        if (A.empty() || B.empty() || A[0].size() != B.size()) {
            return {};
        }
        
        Size m = A.size();
        Size n = B[0].size();
        Size k = A[0].size();
        
        Matrix C(m, Vector(n, 0.0));
        
        for (Size i = 0; i < m; ++i) {
            for (Size j = 0; j < n; ++j) {
                for (Size l = 0; l < k; ++l) {
                    C[i][j] += A[i][l] * B[l][j];
                }
            }
        }
        
        return C;
    }
    
    static Matrix transpose(const Matrix& A) {
        if (A.empty()) return {};
        
        Size m = A.size();
        Size n = A[0].size();
        Matrix AT(n, Vector(m));
        
        for (Size i = 0; i < m; ++i) {
            for (Size j = 0; j < n; ++j) {
                AT[j][i] = A[i][j];
            }
        }
        
        return AT;
    }
    
    static Matrix inverse(const Matrix& A) {
        if (A.empty() || A.size() != A[0].size()) return {};
        
        Size n = A.size();
        Matrix augmented(n, Vector(2 * n, 0.0));
        
        for (Size i = 0; i < n; ++i) {
            for (Size j = 0; j < n; ++j) {
                augmented[i][j] = A[i][j];
            }
            augmented[i][i + n] = 1.0;
        }
        
        for (Size i = 0; i < n; ++i) {
            Size maxRow = i;
            for (Size k = i + 1; k < n; ++k) {
                if (std::abs(augmented[k][i]) > std::abs(augmented[maxRow][i])) {
                    maxRow = k;
                }
            }
            
            if (maxRow != i) {
                std::swap(augmented[i], augmented[maxRow]);
            }
            
            if (std::abs(augmented[i][i]) < EPSILON) return {};
            
            Real pivot = augmented[i][i];
            for (Size j = 0; j < 2 * n; ++j) {
                augmented[i][j] /= pivot;
            }
            
            for (Size k = 0; k < n; ++k) {
                if (k != i) {
                    Real factor = augmented[k][i];
                    for (Size j = 0; j < 2 * n; ++j) {
                        augmented[k][j] -= factor * augmented[i][j];
                    }
                }
            }
        }
        
        Matrix result(n, Vector(n));
        for (Size i = 0; i < n; ++i) {
            for (Size j = 0; j < n; ++j) {
                result[i][j] = augmented[i][j + n];
            }
        }
        
        return result;
    }
    
    static Real determinant(const Matrix& A) {
        if (A.empty() || A.size() != A[0].size()) return 0.0;
        
        Size n = A.size();
        Matrix B = A;
        Real det = 1.0;
        
        for (Size i = 0; i < n; ++i) {
            Size maxRow = i;
            for (Size k = i + 1; k < n; ++k) {
                if (std::abs(B[k][i]) > std::abs(B[maxRow][i])) {
                    maxRow = k;
                }
            }
            
            if (maxRow != i) {
                std::swap(B[i], B[maxRow]);
                det = -det;
            }
            
            if (std::abs(B[i][i]) < EPSILON) return 0.0;
            
            det *= B[i][i];
            
            for (Size k = i + 1; k < n; ++k) {
                Real factor = B[k][i] / B[i][i];
                for (Size j = i; j < n; ++j) {
                    B[k][j] -= factor * B[i][j];
                }
            }
        }
        
        return det;
    }
};

}