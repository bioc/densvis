python_env <- BasiliskEnvironment(
    "densvis",
    pkgname = "densvis",
    packages = c(
        "python=3.12.10",
        "umap-learn=0.5.9.post2",
        "scikit-learn=1.7.0",
        "numba~=0.62",
        "pynndescent=0.5.13",
        "scipy=1.16.0",
        "numpy=~2.3",
        "llvmlite~=0.45"
    )
)
