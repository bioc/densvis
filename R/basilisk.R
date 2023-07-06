python_env <- BasiliskEnvironment(
    "densvis",
    pkgname = "densvis",
    packages = c(
        "umap-learn=0.5.3",
        "scikit-learn=1.3.0",
        "numba=0.57.1",
        "pynndescent=0.5.10",
        "scipy=1.11.1",
        "numpy=1.24.4"
    )
)
