python_env <- BasiliskEnvironment(
    "densvis",
    pkgname = "densvis",
    packages = c(
        "umap-learn=0.5.0",
        "scikit-learn=0.24.0",
        "numba=0.52.0",
        "pynndescent=0.5.1",
        "scipy=1.6.0",
        "numpy=1.19.5"
    )
)
