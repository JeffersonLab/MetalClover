#include "compact_clover.hpp"
using Vec = std::vector<D>;
struct CompactProbe {
  CompactClover c;
  int workers;
  double mass = 0;
  Vec temporary;
  uint64_t products = 0;
  CompactProbe(uint32_t l, int workers, double strength)
      : c(l, 0, workers, strength), workers(workers),
        temporary(c.w.a.p.outputs()) {
    std::vector<D>().swap(c.inverseEven);
    for (auto *v :
         {&c.links32, &c.evenA32, &c.oddA32, &c.inverseEven32, &c.inverseOdd32})
      std::vector<Cx<float>>().swap(*v);
  }
  CD dot(const Vec &x, const Vec &y) {
    std::vector<CD> parts(workers);
    chunks(c.w.a.p.rows * 12, workers,
           [&](int w, uint32_t first, uint32_t last) {
             CD sum = 0;
             for (size_t i = first; i < last; ++i)
               sum += std::conj(CD(x[i].re, x[i].im)) * CD(y[i].re, y[i].im);
             parts[w] = sum;
           });
    CD sum = 0;
    for (auto z : parts)
      sum += z;
    return sum;
  }
  template <class F> void elements(F f) {
    chunks(c.w.a.p.rows * 12, workers, [&](int, uint32_t first, uint32_t last) {
      for (size_t i = first; i < last; ++i)
        f(i);
    });
  }
  void axpy(Vec &x, const Vec &y, CD a) {
    elements([&](size_t i) {
      auto z = a * CD(y[i].re, y[i].im);
      x[i].re += z.real();
      x[i].im += z.imag();
    });
  }
  void scale(Vec &x, double a) {
    elements([&](size_t i) {
      x[i].re *= a;
      x[i].im *= a;
    });
  }
  void Q(const Vec &x, Vec &y) {
    apply_cpu(c.w.a, x.data(), y.data(), workers);
    ++products;
    elements([&](size_t i) {
      double sign = i % 4 < 2 ? 1 : -1;
      y[i] = {sign * (y[i].re + mass * x[i].re),
              sign * (y[i].im + mass * x[i].im)};
    });
  }
  void H(const Vec &x, Vec &y) {
    Q(x, temporary);
    Q(temporary, y);
  }
};
int main(int argc, char **argv) {
  @autoreleasepool {
    try {
      uint32_t l = 32;
      int steps = 150, workers = 8, seed = 12371;
      bool free = false;
      double tol = 1e-5;
      std::vector<double> masses{-.10, -.12, -.14};
      std::string csv = "results/compact-critical32.csv";
      for (int i = 1; i < argc; ++i) {
        std::string f = argv[i];
        auto next = [&]() {
          if (++i >= argc)
            throw std::runtime_error("missing value");
          return std::string(argv[i]);
        };
        if (f == "--size")
          l = std::stoul(next());
        else if (f == "--steps")
          steps = std::stoi(next());
        else if (f == "--seed")
          seed = std::stoi(next());
        else if (f == "--csv")
          csv = next();
        else if (f == "--free")
          free = true;
        else if (f == "--mass") {
          masses.clear();
          std::stringstream s(next());
          std::string x;
          while (std::getline(s, x, ','))
            masses.push_back(std::stod(x));
        } else
          throw std::runtime_error("unknown flag");
      }
      if (l < 4 || l > 32 || l % 2 || steps < 1)
        throw std::runtime_error("invalid probe");
      CompactProbe p(l, workers, free ? 0 : .25);
      size_t n = p.c.w.a.p.outputs();
      std::ofstream out(csv);
      if (!out)
        throw std::runtime_error("cannot open probe CSV");
      out << "extent,mass,sigma_estimate,true_H2_residual,stationary_mass,"
             "steps,converged,products,seconds,seed\n";
      for (double mass : masses) {
        p.mass = mass;
        p.products = 0;
        auto start = Clock::now();
        Vec x(n), ax(n), r(n), ar(n), direction(n), ad(n);
        std::mt19937_64 rng(seed);
        std::normal_distribution<double> rnd(0, 1);
        std::array<CD, 24> constants;
        for (auto &z : constants)
          z = {rnd(rng), rnd(rng)};
        for (size_t i = 0; i < n; ++i) {
          double theta = std::acos(-1.) * double(i / 12 / (l * l * l)) / l;
          auto z = std::cos(theta) * constants[i % 12] +
                   std::sin(theta) * constants[12 + i % 12] +
                   .02 * CD(rnd(rng), rnd(rng));
          x[i] = {z.real(), z.imag()};
        }
        p.scale(x, 1 / std::sqrt(p.dot(x, x).real()));
        p.H(x, ax);
        int it = 0;
        double residual = 0, lambda = 0;
        for (; it < steps; ++it) {
          lambda = p.dot(x, ax).real();
          p.elements([&](size_t i) {
            r[i] = {ax[i].re - lambda * x[i].re, ax[i].im - lambda * x[i].im};
          });
          residual = std::sqrt(p.dot(r, r).real());
          if (it % 25 == 0)
            std::cout << "PROBE m=" << mass << " step=" << it
                      << " sigma=" << std::sqrt(lambda)
                      << " residual=" << residual << std::endl;
          if (residual <= tol)
            break;
          for (int pass = 0; pass < 2; ++pass) {
            auto a = p.dot(x, r);
            p.axpy(r, x, -a);
          }
          double norm = std::sqrt(p.dot(r, r).real());
          if (!(norm > 0))
            throw std::runtime_error("probe residual rank failure");
          p.scale(r, 1 / norm);
          p.H(r, ar);
          for (int pass = 0; pass < 2; ++pass) {
            auto a = p.dot(x, direction);
            p.axpy(direction, x, -a);
            p.axpy(ad, ax, -a);
            auto b = p.dot(r, direction);
            p.axpy(direction, r, -b);
            p.axpy(ad, ar, -b);
          }
          double dn =
              std::sqrt(std::max(0., p.dot(direction, direction).real()));
          int k = dn > 1e-10 ? 3 : 2;
          if (k == 3) {
            p.scale(direction, 1 / dn);
            p.scale(ad, 1 / dn);
          }
          Vec *v[3] = {&x, &r, &direction};
          Vec *av[3] = {&ax, &ar, &ad};
          std::vector<__LAPACK_double_complex> matrix(size_t(k) * k);
          for (int j = 0; j < k; ++j)
            for (int i = 0; i < k; ++i) {
              auto z = p.dot(*v[i], *av[j]);
              matrix[i + k * j] = __LAPACK_double_complex(z);
            }
          __LAPACK_int dim = k, lda = k, info = 0, lwork = 32;
          char job = 'V', uplo = 'U';
          double eigen[3], rw[7];
          __LAPACK_double_complex work[32];
          zheev_(&job, &uplo, &dim, matrix.data(), &lda, eigen, work, &lwork,
                 rw, &info);
          if (info)
            throw std::runtime_error("small probe eigenproblem");
          CD c0 = matrix[0], c1 = matrix[1], c2 = k == 3 ? CD(matrix[2]) : CD{};
          p.elements([&](size_t i) {
            auto d = c1 * CD(r[i].re, r[i].im) +
                     c2 * CD(direction[i].re, direction[i].im),
                 da = c1 * CD(ar[i].re, ar[i].im) + c2 * CD(ad[i].re, ad[i].im);
            auto xx = c0 * CD(x[i].re, x[i].im) + d,
                 aa = c0 * CD(ax[i].re, ax[i].im) + da;
            direction[i] = {d.real(), d.imag()};
            ad[i] = {da.real(), da.imag()};
            x[i] = {xx.real(), xx.imag()};
            ax[i] = {aa.real(), aa.imag()};
          });
          double xn = std::sqrt(p.dot(x, x).real());
          p.scale(x, 1 / xn);
          p.scale(ax, 1 / xn);
        }
        p.H(x, ax);
        lambda = p.dot(x, ax).real();
        p.axpy(ax, x, -lambda);
        residual = std::sqrt(p.dot(ax, ax).real());
        p.Q(x, ax);
        double sigma = std::sqrt(p.dot(ax, ax).real());
        p.elements([&](size_t i) {
          if (i % 4 >= 2) {
            ax[i].re = -ax[i].re;
            ax[i].im = -ax[i].im;
          }
        });
        double stationary = mass - p.dot(x, ax).real();
        if (free) {
          double expected = 1e100;
          for (uint32_t i = 0; i < l * l * l * l; ++i) {
            uint32_t q = i;
            double a = mass, b = 0;
            for (int mu = 0; mu < 4; ++mu) {
              double t = 2 * std::acos(-1.) * (q % l + (mu == 3 ? .5 : 0)) / l;
              q /= l;
              a += 1 - std::cos(t);
              b += std::sin(t) * std::sin(t);
            }
            expected = std::min(expected, std::sqrt(a * a + b));
          }
          if (std::abs(sigma - expected) > 1e-7)
            throw std::runtime_error("free probe analytic gap mismatch");
        }
        if (!std::isfinite(residual) || !std::isfinite(sigma))
          throw std::runtime_error("nonfinite probe");
        out << std::setprecision(14) << l << ',' << mass << ',' << sigma << ','
            << residual << ',' << stationary << ',' << it << ','
            << (residual <= tol) << ',' << p.products << ',' << seconds(start)
            << ',' << seed << '\n';
        out.flush();
        std::cout << "GAP m=" << mass << " estimate=" << sigma
                  << " residual=" << residual << " stationary=" << stationary
                  << " steps=" << it << std::endl;
      }
      return 0;
    } catch (const std::exception &e) {
      std::cerr << "ERROR: " << e.what() << std::endl;
      return 1;
    }
  }
}
