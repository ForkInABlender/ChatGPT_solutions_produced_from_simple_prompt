"""
Single-file SciPy compatibility shim.

Hard constraint:
    This file NEVER imports the real SciPy package.

The scipy.* namespace is constructed dynamically from this file.

Current implemented namespaces:
    scipy.constants
    scipy.special
    scipy.integrate
    scipy.linalg
    scipy.ndimage
    scipy.cluster.vq
    scipy.stats
    scipy.optimize
    scipy.weave

Dependencies:
    - Python standard library
    - NumPy
    - Numba (JIT compilation)

No SciPy dependency exists.
"""

from __future__ import annotations

import math
import sys
import types

import numba
import numba.typed
import numpy as _np

from numpy import (
    __version__ as __numpy_version__,
    reshape, outer, dot, where, array, exp, zeros, size,
    mat, ndarray, eye, tanh, clip, log, sqrt, power, pi, tan, diag,
    random as rand, random, asarray, mgrid, tile, floor, sum,
)
from numpy import *  # noqa: F401,F403  (project-wide star-import kept intentionally)


# ---------------------------------------------------------------------------
# Package identity
# ---------------------------------------------------------------------------

__version__ = "1.1.0-shim"
__all__ = [
    "constants",
    "special",
    "integrate",
    "linalg",
    "ndimage",
    "cluster",
    "stats",
    "optimize",
    "weave",
]


# ---------------------------------------------------------------------------
# Namespace construction  (pure Python — must NOT be JIT-compiled)
# ---------------------------------------------------------------------------

_SHIM_FILE = __file__


def _make_submodule(name: str) -> types.ModuleType:
    """Create a scipy.<name> module whose implementation is this file."""
    fullname = "scipy." + name
    module = types.ModuleType(fullname)
    module.__package__ = "scipy"
    module.__file__ = _SHIM_FILE
    module.__loader__ = globals().get("__loader__")
    module.__spec__ = None
    sys.modules[fullname] = module
    globals()[name] = module
    return module


def _export(module: types.ModuleType, **objects) -> None:
    for name, value in objects.items():
        setattr(module, name, value)


# ---------------------------------------------------------------------------
# scipy.constants
# ---------------------------------------------------------------------------

constants = _make_submodule("constants")

_export(
    constants,
    pi=math.pi,
    e=math.e,
    golden=(1.0 + math.sqrt(5.0)) / 2.0,
    c=299_792_458.0,
    speed_of_light=299_792_458.0,
    h=6.62607015e-34,
    Planck=6.62607015e-34,
    k=1.380649e-23,
    Boltzmann=1.380649e-23,
    N_A=6.02214076e23,
    Avogadro=6.02214076e23,
    R=8.31446261815324,
    gas_constant=8.31446261815324,
    G=6.67430e-11,
    gravitational_constant=6.67430e-11,
    epsilon_0=8.8541878128e-12,
    mu_0=1.25663706212e-6,
    elementary_charge=1.602176634e-19,
    eV=1.602176634e-19,
    electron_mass=9.1093837015e-31,
    proton_mass=1.67262192369e-27,
    neutron_mass=1.67492749804e-27,
    alpha=7.2973525693e-3,
)


# ---------------------------------------------------------------------------
# scipy.special
# ---------------------------------------------------------------------------

special = _make_submodule("special")

# ------------------------------------------------------------------
# Scalar JIT helpers (nopython=True, cache=True)
# ------------------------------------------------------------------

@numba.jit(nopython=True, cache=True)
def _gamma(x):
    return math.gamma(x)


@numba.jit(nopython=True, cache=True)
def _gammaln(x):
    return math.lgamma(x)


@numba.jit(nopython=True, cache=True)
def _erf(x):
    return math.erf(x)


@numba.jit(nopython=True, cache=True)
def _erfc(x):
    return math.erfc(x)


@numba.jit(nopython=True, cache=True)
def _expit(x):
    return 1.0 / (1.0 + math.exp(-x))


@numba.jit(nopython=True, cache=True)
def _log_expit(x):
    if x >= 0:
        return -math.log1p(math.exp(-x))
    return x - math.log1p(math.exp(x))


@numba.jit(nopython=True, cache=True)
def _logit(x):
    if x <= 0.0 or x >= 1.0:
        return math.nan
    return math.log(x / (1.0 - x))


@numba.jit(nopython=True, cache=True)
def _logsumexp(values):
    """
    Log-sum-exp over a 1-D float64 array (nopython, no boxing).
    """
    n = len(values)
    if n == 0:
        return -math.inf

    m = values[0]
    for i in range(1, n):
        if values[i] > m:
            m = values[i]

    if math.isinf(m):
        return m

    acc = 0.0
    for i in range(n):
        acc += math.exp(values[i] - m)

    return m + math.log(acc)


@numba.jit(nopython=True, cache=True)
def _comb(n, k):
    return math.comb(n, k)


@numba.jit(nopython=True, cache=True)
def _perm(n, k):
    return math.perm(n, k)


@numba.jit(nopython=True, cache=True)
def _factorial(n):
    return math.factorial(n)


# ------------------------------------------------------------------
# Vectorized ufuncs — compile to native ufunc over float64 arrays.
# These accept scalars *and* arrays without Python overhead per element.
# ------------------------------------------------------------------

@numba.vectorize(["float64(float64)"], nopython=True, cache=True)
def _vec_erf(x):
    return math.erf(x)


@numba.vectorize(["float64(float64)"], nopython=True, cache=True)
def _vec_erfc(x):
    return math.erfc(x)


@numba.vectorize(["float64(float64)"], nopython=True, cache=True)
def _vec_gamma(x):
    return math.gamma(x)


@numba.vectorize(["float64(float64)"], nopython=True, cache=True)
def _vec_gammaln(x):
    return math.lgamma(x)


@numba.vectorize(["float64(float64)"], nopython=True, cache=True)
def _vec_expit(x):
    return 1.0 / (1.0 + math.exp(-x))


@numba.vectorize(["float64(float64)"], nopython=True, cache=True)
def _vec_logit(x):
    if x <= 0.0 or x >= 1.0:
        return math.nan
    return math.log(x / (1.0 - x))


@numba.vectorize(["float64(float64)"], nopython=True, cache=True)
def _vec_log_expit(x):
    if x >= 0:
        return -math.log1p(math.exp(-x))
    return x - math.log1p(math.exp(x))


# Public wrappers: scalar path calls the JIT scalar; array path calls vectorize.

def gamma(x):
    if _np.ndim(x) == 0:
        return _gamma(float(x))
    return _vec_gamma(_np.asarray(x, dtype=float))


def gammaln(x):
    if _np.ndim(x) == 0:
        return _gammaln(float(x))
    return _vec_gammaln(_np.asarray(x, dtype=float))


def erf(x):
    if _np.ndim(x) == 0:
        return _erf(float(x))
    return _vec_erf(_np.asarray(x, dtype=float))


def erfc(x):
    if _np.ndim(x) == 0:
        return _erfc(float(x))
    return _vec_erfc(_np.asarray(x, dtype=float))


def expit(x):
    if _np.ndim(x) == 0:
        return _expit(float(x))
    return _vec_expit(_np.asarray(x, dtype=float))


def log_expit(x):
    if _np.ndim(x) == 0:
        return _log_expit(float(x))
    return _vec_log_expit(_np.asarray(x, dtype=float))


def logit(x):
    if _np.ndim(x) == 0:
        return _logit(float(x))
    return _vec_logit(_np.asarray(x, dtype=float))


def logsumexp(a, axis=None, b=None, keepdims=False, return_sign=False):
    """
    Log-sum-exp with optional axis/weights, matching scipy.special.logsumexp.
    1-D no-axis path uses the fast Numba kernel; all other paths use NumPy.
    """
    a = _np.asarray(a, dtype=float)

    if b is not None:
        b = _np.asarray(b, dtype=float)

    if axis is None and not keepdims and b is None and not return_sign:
        # Fast path: flatten and use the nopython kernel.
        result = _logsumexp(a.ravel())
        if return_sign:
            return result, 1.0
        return result

    # General NumPy path.
    if b is not None:
        a_max = _np.max(a, axis=axis, keepdims=True)
        tmp = b * _np.exp(a - a_max)
        s = _np.sum(tmp, axis=axis, keepdims=keepdims)
        if return_sign:
            sgn = _np.sign(s)
            out = _np.log(_np.abs(s))
            out += _np.squeeze(a_max, axis=axis) if not keepdims else a_max
            return out, sgn
        out = _np.log(s)
        out += _np.squeeze(a_max, axis=axis) if not keepdims else a_max
        return out

    a_max = _np.max(a, axis=axis, keepdims=True)
    tmp = _np.exp(a - a_max)
    s = _np.sum(tmp, axis=axis, keepdims=keepdims)
    out = _np.log(s)
    out += _np.squeeze(a_max, axis=axis) if not keepdims else a_max
    if return_sign:
        return out, _np.ones_like(out)
    return out


def comb(n, k):
    return _comb(int(n), int(k))


def perm(n, k=None):
    if k is None:
        k = n
    return _perm(int(n), int(k))


def factorial(n):
    return _factorial(int(n))


_export(
    special,
    gamma=gamma,
    gammaln=gammaln,
    erf=erf,
    erfc=erfc,
    expit=expit,
    log_expit=log_expit,
    logit=logit,
    logsumexp=logsumexp,
    comb=comb,
    perm=perm,
    factorial=factorial,
)


# ---------------------------------------------------------------------------
# scipy.integrate
# ---------------------------------------------------------------------------

integrate = _make_submodule("integrate")


def quad(
    func,
    a,
    b,
    args=(),
    epsabs=1.49e-8,
    epsrel=1.49e-8,
    limit=50,
):
    """
    Numerically integrate func from a to b using adaptive composite Simpson.

    Returns:
        (integral, estimated_error)
    """
    if args is None:
        args = ()

    if limit < 1:
        raise ValueError("limit must be positive")

    n = 32
    previous = None
    result = None

    while True:
        if n % 2:
            n += 1

        h = (b - a) / n
        total = func(a, *args) + func(b, *args)

        for i in range(1, n):
            x = a + i * h
            if i & 1:
                total += 4.0 * func(x, *args)
            else:
                total += 2.0 * func(x, *args)

        result = total * h / 3.0

        if previous is not None:
            error = abs(result - previous) / 15.0
            tolerance = max(epsabs, epsrel * abs(result))
            if error <= tolerance:
                return result, error

        if n >= max(2, limit * 2):
            if previous is None:
                return result, 0.0
            return result, abs(result - previous) / 15.0

        previous = result
        n *= 2


def simpson(y, x=None, dx=1.0, axis=-1):
    """
    Composite Simpson integration for sampled data.
    """
    y = _np.asarray(y)

    if y.shape[axis] < 2:
        raise ValueError("at least two samples are required")

    y = _np.moveaxis(y, axis, -1)
    n = y.shape[-1]

    if x is None:
        spacing = dx
    else:
        x = _np.asarray(x)
        if x.ndim != 1 or x.size != n:
            raise ValueError(
                "x must have the same length as the integration axis"
            )
        spacing = None

    if n == 2:
        if x is None:
            return dx * (y[..., 0] + y[..., 1]) / 2.0
        return (x[1] - x[0]) * (y[..., 0] + y[..., 1]) / 2.0

    if n % 2 == 1:
        if x is None:
            h = spacing
            return (
                y[..., 0]
                + y[..., -1]
                + 4.0 * _np.sum(y[..., 1:-1:2], axis=-1)
                + 2.0 * _np.sum(y[..., 2:-2:2], axis=-1)
            ) * h / 3.0

        h = _np.diff(x)
        total = _np.zeros(y.shape[:-1], dtype=_np.result_type(y, float))
        for i in range(0, n - 2, 2):
            h0 = h[i]
            h1 = h[i + 1]
            total += (
                (h0 + h1) / 6.0
                * (
                    (2.0 - h1 / h0) * y[..., i]
                    + ((h0 + h1) ** 2 / (h0 * h1)) * y[..., i + 1]
                    + (2.0 - h0 / h1) * y[..., i + 2]
                )
            )
        return total

    # Even number of samples: Simpson over all but last interval, trapezoid for last.
    if x is None:
        h = spacing
        result = (
            y[..., 0]
            + y[..., -2]
            + 4.0 * _np.sum(y[..., 1:-2:2], axis=-1)
            + 2.0 * _np.sum(y[..., 2:-3:2], axis=-1)
        ) * h / 3.0
        result += h * (y[..., -2] + y[..., -1]) / 2.0
        return result

    return simpson(y[..., :-1], x=x[:-1], axis=-1) + (
        x[-1] - x[-2]
    ) * (y[..., -2] + y[..., -1]) / 2.0


_export(integrate, quad=quad, simpson=simpson)


# ---------------------------------------------------------------------------
# scipy.linalg
# ---------------------------------------------------------------------------

linalg = _make_submodule("linalg")


@numba.jit(nopython=True, cache=True)
def _solve_nb(a, b):
    """Thin nopython wrapper so linalg.solve benefits from Numba dispatch."""
    return _np.linalg.solve(a, b)


def solve(a, b, **kwargs):  # **kwargs for scipy compatibility
    return _solve_nb(_np.asarray(a, dtype=float), _np.asarray(b, dtype=float))


def inv(a):
    return _np.linalg.inv(a)


def det(a):
    return _np.linalg.det(a)


def slogdet(a):
    return _np.linalg.slogdet(a)


def eig(a):
    return _np.linalg.eig(a)


def eigh(a, UPLO="L"):
    return _np.linalg.eigh(a, UPLO=UPLO)


def eigvals(a):
    return _np.linalg.eigvals(a)


def eigvalsh(a, UPLO="L"):
    return _np.linalg.eigvalsh(a, UPLO=UPLO)


def svd(a, full_matrices=True, compute_uv=True, hermitian=False):
    return _np.linalg.svd(
        a,
        full_matrices=full_matrices,
        compute_uv=compute_uv,
        hermitian=hermitian,
    )


def norm(a, ord=None, axis=None, keepdims=False):
    return _np.linalg.norm(a, ord=ord, axis=axis, keepdims=keepdims)


def matrix_power(a, n):
    return _np.linalg.matrix_power(a, n)


def matrix_rank(a, tol=None, hermitian=False):
    return _np.linalg.matrix_rank(a, tol=tol, hermitian=hermitian)


def pinv(a, **kwargs):
    return _np.linalg.pinv(a)


# pinv2 is scipy's legacy alias
pinv2 = pinv


def lstsq(a, b, rcond=None):
    return _np.linalg.lstsq(a, b, rcond=rcond)


def cholesky(a, upper=False):
    result = _np.linalg.cholesky(a)
    return result.T.conj() if upper else result


def qr(a, mode="reduced"):
    return _np.linalg.qr(a, mode=mode)


def orth(A):
    """Return an orthonormal basis for the range of A."""
    Q, R = _np.linalg.qr(A, mode="reduced")
    tol = max(A.shape) * _np.finfo(float).eps * abs(R).max()
    rank = _np.sum(_np.abs(_np.diag(R)) > tol)
    return Q[:, :rank]


_export(
    linalg,
    solve=solve,
    inv=inv,
    det=det,
    slogdet=slogdet,
    eig=eig,
    eigh=eigh,
    eigvals=eigvals,
    eigvalsh=eigvalsh,
    svd=svd,
    norm=norm,
    matrix_power=matrix_power,
    matrix_rank=matrix_rank,
    pinv=pinv,
    pinv2=pinv2,
    lstsq=lstsq,
    cholesky=cholesky,
    qr=qr,
    orth=orth,
)


# ---------------------------------------------------------------------------
# scipy.ndimage
# ---------------------------------------------------------------------------

ndimage = _make_submodule("ndimage")


def minimum_position(input):
    """Return the (row, col, …) position of the minimum value."""
    inp = _np.asarray(input)
    idx = _np.argmin(inp)
    return tuple(int(i) for i in _np.unravel_index(idx, inp.shape))


_export(ndimage, minimum_position=minimum_position)


# ---------------------------------------------------------------------------
# scipy.cluster.vq
# ---------------------------------------------------------------------------

cluster = _make_submodule("cluster")
cluster_vq = types.ModuleType("scipy.cluster.vq")
cluster_vq.__package__ = "scipy.cluster"
cluster_vq.__file__ = _SHIM_FILE
cluster_vq.__spec__ = None
sys.modules["scipy.cluster.vq"] = cluster_vq
setattr(cluster, "vq", cluster_vq)


@numba.jit(nopython=True, cache=True, parallel=True)
def _kmeans_assign(data, centroids):
    """
    Parallel assignment step: returns label array (int64).
    data      : (n, d) float64
    centroids : (k, d) float64
    """
    n = data.shape[0]
    k = centroids.shape[0]
    labels = _np.empty(n, dtype=_np.int64)
    for i in numba.prange(n):
        best_j = 0
        best_d = 0.0
        for dim in range(data.shape[1]):
            diff = data[i, dim] - centroids[0, dim]
            best_d += diff * diff
        for j in range(1, k):
            d = 0.0
            for dim in range(data.shape[1]):
                diff = data[i, dim] - centroids[j, dim]
                d += diff * diff
            if d < best_d:
                best_d = d
                best_j = j
        labels[i] = best_j
    return labels


@numba.jit(nopython=True, cache=True)
def _kmeans_update(data, labels, k):
    """
    Update centroids as the mean of assigned points.
    Returns new centroids (k, d) float64.
    """
    d = data.shape[1]
    centroids = _np.zeros((k, d))
    counts = _np.zeros(k, dtype=_np.int64)
    for i in range(data.shape[0]):
        lbl = labels[i]
        counts[lbl] += 1
        for dim in range(d):
            centroids[lbl, dim] += data[i, dim]
    for j in range(k):
        if counts[j] > 0:
            for dim in range(d):
                centroids[j, dim] /= counts[j]
    return centroids


def kmeans2(data, k, iter=10, minit="random", **kwargs):
    """
    Pure-NumPy/Numba k-means clustering (subset of scipy.cluster.vq.kmeans2).

    Returns:
        centroid : (k, d) float64 array of cluster centres
        label    : (n,)   int64 array of assignments
    """
    data = _np.asarray(data, dtype=_np.float64)
    scalar_input = data.ndim == 1
    if scalar_input:
        data = data[:, None]

    n, d = data.shape
    rng = _np.random.default_rng()

    if isinstance(k, int):
        idx = rng.choice(n, size=k, replace=False)
        centroids = data[idx].copy()
    else:
        centroids = _np.asarray(k, dtype=_np.float64)
        k = len(centroids)

    labels = _np.zeros(n, dtype=_np.int64)

    for _ in range(iter):
        labels = _kmeans_assign(data, centroids)
        new_centroids = _kmeans_update(data, labels, k)
        if _np.allclose(centroids, new_centroids):
            break
        centroids = new_centroids

    if scalar_input:
        centroids = centroids.ravel()

    return centroids, labels


_export(cluster_vq, kmeans2=kmeans2)


# ---------------------------------------------------------------------------
# scipy.stats
# ---------------------------------------------------------------------------

stats = _make_submodule("stats")

# --- Numba-accelerated normal PPF (Beasley-Springer-Moro) ---

@numba.vectorize(["float64(float64)"], nopython=True, cache=True)
def _norm_ppf_scalar(q):
    """Per-element inverse-normal CDF (Beasley-Springer-Moro)."""
    if q <= 0.0:
        return -math.inf
    if q >= 1.0:
        return math.inf
    p = q if q < 0.5 else 1.0 - q
    t = math.sqrt(-2.0 * math.log(p))
    num = 2.515517 + 0.802853 * t + 0.010328 * t * t
    den = 1.0 + 1.432788 * t + 0.189269 * t * t + 0.001308 * t * t * t
    z = t - num / den
    return -z if q < 0.5 else z


@numba.vectorize(["float64(float64)"], nopython=True, cache=True)
def _norm_cdf_scalar(x):
    """Per-element standard normal CDF via erf."""
    return 0.5 * (1.0 + math.erf(x * 0.7071067811865476))  # 1/sqrt(2)


class _NormDist:
    """
    Minimal frozen/un-frozen normal distribution matching scipy.stats.norm.
    Vectorised paths use Numba ufuncs; no Python loops over array elements.
    """

    def pdf(self, x, loc=0.0, scale=1.0):
        x = _np.asarray(x, dtype=float)
        return _np.exp(-0.5 * ((x - loc) / scale) ** 2) / (
            scale * _np.sqrt(2 * _np.pi)
        )

    def logpdf(self, x, loc=0.0, scale=1.0):
        x = _np.asarray(x, dtype=float)
        return -0.5 * ((x - loc) / scale) ** 2 - _np.log(
            scale * _np.sqrt(2 * _np.pi)
        )

    def cdf(self, x, loc=0.0, scale=1.0):
        """CDF via Numba vectorize — no Python loop, no np.vectorize."""
        x = _np.asarray(x, dtype=float)
        return _norm_cdf_scalar((x - loc) / scale)

    def ppf(self, q, loc=0.0, scale=1.0):
        """Percent-point via Numba vectorize ufunc — O(1) Python overhead."""
        q = _np.asarray(q, dtype=float)
        scalar = q.ndim == 0
        result = _norm_ppf_scalar(_np.atleast_1d(q))
        result = loc + scale * result
        return float(result[0]) if scalar else result

    def rvs(self, loc=0.0, scale=1.0, size=None):
        return _np.random.normal(loc=loc, scale=scale, size=size)

    def __call__(self, loc=0.0, scale=1.0):
        """Return a frozen instance with fixed loc/scale."""
        parent = self

        class _Frozen:
            def pdf(self, x):
                return parent.pdf(x, loc, scale)

            def logpdf(self, x):
                return parent.logpdf(x, loc, scale)

            def cdf(self, x):
                return parent.cdf(x, loc, scale)

            def ppf(self, q):
                return parent.ppf(q, loc, scale)

            def rvs(self, size=None):
                return parent.rvs(loc, scale, size)

        return _Frozen()


_export(stats, norm=_NormDist())


# ---------------------------------------------------------------------------
# scipy.optimize
# ---------------------------------------------------------------------------

optimize = _make_submodule("optimize")


def fmin(
    func,
    x0,
    args=(),
    xtol=1e-4,
    ftol=1e-4,
    maxiter=None,
    maxfun=None,
    full_output=False,
    disp=True,
    retall=False,
    callback=None,
    **kwargs,
):
    """
    Nelder-Mead simplex minimization (mirrors scipy.optimize.fmin).
    """
    x0 = _np.asarray(x0, dtype=float).ravel()
    n = len(x0)

    if maxiter is None:
        maxiter = n * 200
    if maxfun is None:
        maxfun = n * 200

    # Build initial simplex
    sim = _np.zeros((n + 1, n))
    sim[0] = x0.copy()
    for i in range(n):
        x = x0.copy()
        x[i] += 0.05 if x[i] != 0 else 0.00025
        sim[i + 1] = x

    fvals = _np.array([func(sim[i], *args) for i in range(n + 1)])
    fcalls = n + 1
    iterations = 0
    allvecs = [sim[0].copy()] if retall else None

    rho, chi, psi, sigma = 1.0, 2.0, 0.5, 0.5

    while fcalls < maxfun and iterations < maxiter:
        order = _np.argsort(fvals)
        sim = sim[order]
        fvals = fvals[order]

        if retall:
            allvecs.append(sim[0].copy())

        # Convergence check
        if max(abs(fvals[1:] - fvals[0])) <= ftol and max(
            _np.max(abs(sim[1:] - sim[0]), axis=1)
        ) <= xtol:
            break

        xbar = sim[:-1].mean(axis=0)
        xr = (1 + rho) * xbar - rho * sim[-1]
        fr = func(xr, *args)
        fcalls += 1

        if fr < fvals[0]:
            xe = (1 + rho * chi) * xbar - rho * chi * sim[-1]
            fe = func(xe, *args)
            fcalls += 1
            if fe < fr:
                sim[-1] = xe
                fvals[-1] = fe
            else:
                sim[-1] = xr
                fvals[-1] = fr
        elif fr < fvals[-2]:
            sim[-1] = xr
            fvals[-1] = fr
        else:
            if fr < fvals[-1]:
                xc = (1 + psi * rho) * xbar - psi * rho * sim[-1]
                fc = func(xc, *args)
                fcalls += 1
                if fc <= fr:
                    sim[-1] = xc
                    fvals[-1] = fc
                else:
                    sim[1:] = sim[0] + sigma * (sim[1:] - sim[0])
                    fvals[1:] = [func(sim[i], *args) for i in range(1, n + 1)]
                    fcalls += n
            else:
                xc = (1 - psi) * xbar + psi * sim[-1]
                fc = func(xc, *args)
                fcalls += 1
                if fc < fvals[-1]:
                    sim[-1] = xc
                    fvals[-1] = fc
                else:
                    sim[1:] = sim[0] + sigma * (sim[1:] - sim[0])
                    fvals[1:] = [func(sim[i], *args) for i in range(1, n + 1)]
                    fcalls += n

        if callback is not None:
            callback(sim[0])

        iterations += 1

    x = sim[0]
    fval = fvals[0]

    if full_output:
        if retall:
            return x, fval, iterations, fcalls, 0, allvecs
        return x, fval, iterations, fcalls, 0
    if retall:
        return x, allvecs
    return x


_export(optimize, fmin=fmin)


# ---------------------------------------------------------------------------
# scipy.weave  (removed from scipy; stub to avoid ImportError)
# ---------------------------------------------------------------------------

weave = _make_submodule("weave")


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------

def _verify_no_external_scipy():
    """
    Verify that every scipy.* module currently loaded belongs to this shim.
    """
    for name, module in list(sys.modules.items()):
        if not name.startswith("scipy"):
            continue
        if name != "scipy" and not name.startswith("scipy."):
            raise RuntimeError("Unexpected module in scipy namespace: " + name)
        module_file = getattr(module, "__file__", None)
        if module_file is not None and module_file != _SHIM_FILE:
            raise RuntimeError(
                "External SciPy module detected: " + name + " -> " + str(module_file)
            )


# _verify_no_external_scipy()
