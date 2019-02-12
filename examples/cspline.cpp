// clang++ cspline.cpp -o cspline --std=c++2a -stdlib=libc++

#include <cstdio>
#include <limits>
#include <optional>
#include <vector>
#include <string>
#include <fstream>

#define DEBUG_VERBOSE 0

using number_t = float;
const number_t NaN = std::numeric_limits<number_t>::quiet_NaN();

// p: (x, y)
struct Point2 {
    number_t x;
    number_t y;
};

// S: a(x-xs)^3 + b(x-xs)^2 + c(x-xs) + d
struct SectionCurve2 {
    number_t a;
    number_t b;
    number_t c;
    number_t d;
    number_t xs;
    number_t xe;
};

static std::vector<Point2> read_points_from_csv_file(char const * path);
static std::vector<SectionCurve2> cspline_for_points(std::vector<Point2> const & points, std::optional<number_t> opt_a0);
static number_t y_on_cspline(std::vector<SectionCurve2> const & curves, number_t x);
static number_t y_on_curve(SectionCurve2 const & curve, number_t x);

int main(int argc, char ** argv) {
    if (argc < 2) {
        std::puts("usage: cspline source [a0]");
        return 0;
    }
    std::optional<number_t> opt_a0;
    if (argc >= 3) {
        number_t a0;
        if (std::sscanf(argv[2], "%f", &a0) != 1) {
            std::printf("Error %d\n", -__LINE__);
            return -1;
        }
        opt_a0.emplace(a0);
    }
    try {
        // input
        auto const points = read_points_from_csv_file(argv[1]);
        auto const curves = cspline_for_points(points, opt_a0);
        // output
        {
            auto range = points.back().x - points.front().x;
            for (auto x = points.front().x; x < points.back().x; x += range / 100) {
                std::printf("%f,%f\n", x, y_on_cspline(curves, x));
            }
            std::printf("%f,%f\n", points.back().x, y_on_cspline(curves, points.back().x));
        }
    } catch (int & error) {
        std::printf("Error %d\n", error);
        return -1;
    }

    return 0;
}

static std::vector<Point2> read_points_from_csv_file(char const * path) {
    std::vector<Point2> points;

    std::ifstream ifs(path);
    if (ifs.fail()) {
        throw -__LINE__;
    }

    std::string line;
    while (std::getline(ifs, line)) {
        decltype(points)::value_type point;
        if (std::sscanf(line.c_str(), "%f,%f", &point.x, &point.y) != 2) {
            throw -__LINE__;
        }
        points.push_back(point);
    }
    return points;
}

static std::vector<SectionCurve2> cspline_for_points(std::vector<Point2> const & points, std::optional<number_t> opt_a0) {
    std::vector<SectionCurve2> S;

    if (!(points.size() >= 2)) {
        throw -__LINE__;
    }
    for (auto i = 0; i < points.size() - 1; ++i) {
        if (!(points[i].x < points[i+1].x)) {
            throw -__LINE__;
        }
    }

    auto const N = points.size() - 1; // N+1 points => N curves

    // Prepare curves
    for (auto i = 0; i < N; ++i) {
        SectionCurve2 curve = {
            .a = NaN,
            .b = NaN,
            .c = NaN,
            .d = NaN,
            .xs = points[i].x,
            .xe = points[i+1].x
        };
        S.push_back(curve);
    }
    // d_i
    for (auto i = 0; i < N; ++i) {
        S[i].d = points[i].y;
    }
    // b_i
    {
        // X
        auto X = new number_t[(N-1) * (N-1)];
        {
            for (auto r = 0; r < N-1; ++r) {
                for (auto c = 0; c < N-1; ++c) {
                    X[r * (N-1) + c] = 0;
                }
            }
            // 0
            auto xd0 = points[1].x - points[0].x;
            auto xd1 = points[2].x - points[1].x;
            X[0 * (N-1) + 0] = 2 * (xd0 + xd1);
            X[0 * (N-1) + 1] = xd1;
            // 1, 2, ... n-3
            for (auto r = 1; r < N-2; ++r) {
                auto c = r-1;
                X[r * (N-1) + c+0] = points[r+1].x - points[r].x;
                X[r * (N-1) + c+1] = 2 * ((points[r+1].x - points[r].x) + (points[r+2].x - points[r+1].x));
                X[r * (N-1) + c+2] = points[r+2].x - points[r+1].x;
            }
            // n-2
            auto xdn_2 = points[N-1].x - points[N-2].x;
            auto xdn_1 = points[N].x - points[N-1].x;
            X[(N-2) * (N-1) + N-3] = xdn_2;
            X[(N-2) * (N-1) + N-2] = 2 * (xdn_2 + xdn_1);

            if (DEBUG_VERBOSE) {
                printf("X = | ");
                for (auto r = 0; r < N-1; ++r) {
                    for (auto c = 0; c < N-1; ++c) {
                        if (c > 0) {
                            printf(" ");
                        }
                        printf("%10.3e", X[r * (N-1) + c]);
                    }
                    if (r != N-2) {
                        printf(" |\n    | ");
                    } else {
                        printf(" |\n");
                    }
                }
            }
        }
        // X^{-1}
        auto X_inv = new number_t[(N-1) * (N-1)];
        {
            for (auto r = 0; r < N-1; ++r) {
                for (auto c = 0; c < N-1; ++c) {
                    X_inv[r * (N-1) + c] = r == c ? 1 : 0;
                }
            }
            // row reduction: bottom-left (destructs X)
            for (auto r = 1; r < N-1; ++r) {
                auto scale = X[r * (N-1) + r-1] / X[(r-1) * (N-1) + r-1];
                for (auto c = 0; c < N-1; ++c) {
                    X[r * (N-1) + c] -= scale * X[(r-1) * (N-1) + c];
                    X_inv[r * (N-1) + c] -= scale * X_inv[(r-1) * (N-1) + c];
                }
            }
            // row reduction: top-right (destructs X)
            for (auto r_inv = 1; r_inv < N-1; ++r_inv) {
                auto r = N-1 - r_inv;
                auto scale = X[(r-1) * (N-1) + r] / X[r * (N-1) + r];
                for (auto c = 0; c < N-1; ++c) {
                    X[(r-1) * (N-1) + c] -= scale * X[r * (N-1) + c];
                    X_inv[(r-1) * (N-1) + c] -= scale * X_inv[r * (N-1) + c];
                }
            }
            // row reduction: normalize
            for (auto r = 0; r < N-1; ++r) {
                auto scale = 1 / X[r * (N-1) + r];
                for (auto c = 0; c < N-1; ++c) {
                    X[r * (N-1) + c] *= scale;
                    X_inv[r * (N-1) + c] *= scale;
                }
            }

            if (DEBUG_VERBOSE) {
                printf("X = | ");
                for (auto r = 0; r < N-1; ++r) {
                    for (auto c = 0; c < N-1; ++c) {
                        if (c > 0) {
                            printf(" ");
                        }
                        printf("%10.3e", X[r * (N-1) + c]);
                    }
                    if (r != N-2) {
                        printf(" |\n    | ");
                    } else {
                        printf(" |\n");
                    }
                }

                printf("X^{-1} = | ");
                for (auto r = 0; r < N-1; ++r) {
                    for (auto c = 0; c < N-1; ++c) {
                        if (c > 0) {
                            printf(" ");
                        }
                        printf("%10.3e", X_inv[r * (N-1) + c]);
                    }
                    if (r != N-2) {
                        printf(" |\n         | ");
                    } else {
                        printf(" |\n");
                    }
                }
            }
        }
        auto Z = new number_t[N-1];
        for (auto i = 0; i < N-1; ++i) {
            auto ydi1 = points[i+2].y - points[i+1].y;
            auto xdi1 = points[i+2].x - points[i+1].x;
            auto ydi = points[i+1].y - points[i].y;
            auto xdi = points[i+1].x - points[i].x;
            Z[i] = 3 * (ydi1 / xdi1 - ydi / xdi);
        }
        // b_i
        {
            S[0].b = 0;
            for (auto i = 1; i < N; ++i) {
                S[i].b = 0;
                for (auto j = 0; j < N-1; ++j) {
                    S[i].b += X_inv[(i-1) * (N-1) + j] * Z[j];
                }
            }
        }
        delete[] X;
        delete[] X_inv;
        delete[] Z;
    }
    // a_i
    {
        for (auto i = 0; i < N-1; ++i) {
            auto xdi = points[i+1].x - points[i].x;
            S[i].a = (S[i+1].b - S[i].b) / (3 * xdi);
        }
        static auto const bn = 0;
        auto xdn_1 = points[N].x - points[N-1].x;
        S[N-1].a = (bn - S[N-1].b) / (3 * xdn_1);
    }
    // c_i
    for (auto i = 0; i < N; ++i) {
        auto dy = points[i+1].y - points[i].y;
        auto dx = points[i+1].x - points[i].x;

        S[i].c = dy / dx - S[i].a * dx * dx - S[i].b * dx;
    }
    return S;
}

static number_t y_on_cspline(std::vector<SectionCurve2> const & curves, number_t x) {
    for (auto curve: curves) {
        if (x >= curve.xs && x <= curve.xe) {
            return y_on_curve(curve, x);
        }
    }
    return NaN;
}

static number_t y_on_curve(SectionCurve2 const & curve, number_t x) {
    auto nx = x - curve.xs;
    return curve.a * nx * nx * nx + curve.b * nx * nx + curve.c * nx + curve.d;
}
